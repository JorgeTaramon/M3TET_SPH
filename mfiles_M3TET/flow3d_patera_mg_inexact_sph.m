function [VAR,MESH,Uinfo] = flow3d_patera_mg_inexact(VAR,MESH,UBC,SETTINGS,PHYSICS,NUMSCALE,dt)
% Usage: [VAR,MESH,Uinfo] = flow3d_patera_mg_inexact(VAR,MESH,UBC,SETTINGS,PHYSICS,NUMSCALE,dt)
%
% Purpose: Stokes flow solver (calculate velocity and pressure solution)
%
% Input:
%   VAR      : [structure] : major variable fields, each is a vector
%   MESH     : [structure] : FE mesh parameters
%   UBC      : [structure] : velocity boundary conditions
%   SETTINGS : [structure] : model parameters
%   PHYSICS  : [structure] : physical properties
%   NUMSCALE : [structure] : numerical scaling parameters
%   dt       : [double]    : time step
%
% Output:
%   VAR      : [structure] : major variable fields, each is a vector
%   MESH     : [structure] : FE mesh parameters
%   Uinfo    : [structure] : profiling data
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% Details on the Stokes flow solver:
% - Patera pressure iterations (momentum eq. formally solved for velocity
%   and substituted into incompressibility constraint)
% - non-parallel version
% - elements: 3D Taylor-Hood triangles (10 velocity nodes, 4 pressure nodes)
% - velocity solution: Conjugate Gradient algorithm, preconditioned by a
%                      single V-cylce of geometric multigrid
% - pressure solution: Preconditioned Conjugate Gradient (PCG) solver
%   (for available preconditioners see "Setup"-block below)
% - Null-space removed by subtracting the mean of the pressure solution
%   during the iterations (only if boundary conditions do not allow for 
%   unconstrained in- or outflow)
% - flexible CG used to orthogonalize q allows using inexact (thus slightly
%   asymmetric) preconditioner
% - element assembly using standard method OR much faster blockwise 
%   assembly (MILAMIN style)
% - purely viscous rheology (no elasticity yet)
%
% JH Jan 2011
% JH/JPM Feb2011
% JH Mar/Apr 2011
% JH May 2014: - max. velocity change during iterations now criterium for
%                ending iterations
%              - added free surface (algorithm by Miguel & Jason)
                                                                            clock0=clock;tic

% =========================================================================
% SETUP FOR STOKES FLOW SOLVER
% =========================================================================
FigNo    = 31;
    % number of figure showing convergence; ==0 --> no plot
fs       = 10;
    % font size in convergence plot
YLimit   = [-10 10];
    % y-axis limits in convergence plot

% PART 1 - PRESSURE ITERATIONS
rtol_Pat   = 1e-3;
    % RELATIVE tolerance for pressure solution
    % (relative to norm of initial pressure residual)
PC_Pat     = 4;
    % 0 ==> no preconditioner (d=r)
    % 1 ==> d = r./diag(MM);  scaling by diagonal of mass matrix
    % 2 ==> d = r./LMM; scaling by lumped mass matrix; LMM=sum(MM,2)
    % 3 ==> d = d + (r - MM*d)./LMM; Jacobi iterations on MM
    % 4 ==> inexact Patera algorithm solving S d = r
nit_MM     = 10;
    % number of mass matrix iterations when using PC_Pat==3 or PC_Pat==4
if PC_Pat==4
    itmin_Pat = 1;   % minimum number of pressure iterations
    itmax_Pat = 20;  % maximum number of pressure iterations
else
    itmin_Pat = 3;   % minimum number of pressure iterations
    itmax_Pat = 120; % maximum number of pressure iterations
end

% PART 2 - PRECONDITIONING OF INEXACT PATERA ALGORITHM (PC_Pat==4)
OPTS_PC.rtol   = 2e-2;
    % tolerance of preconditioning (inexact) Patera relative to current 
    % norm of pressure residual
OPTS_PC.itmin  = 3;
    % minimum number of inexact Patera iterations
OPTS_PC.itmax  = 100;
    % maximum number of inexact Patera iterations
OPTS_PC.PC     = 3;
    % preconditioner used by inexact Patera (see list above)
OPTS_PC.rtol_Z = 1e-2;
    % tolerance for K z = y solution inside inexact Patera
    % (relative to current pressure tolerance inside inexact Patera)

% PART 3 - VELOCITY SOLUTION
rtol_U     = 1e-2;
    % tolerance for velocity solution (RELATIVE to pressure tolerance; use ~1e-2)
itmax_U    = 30;
    % maximum number of velocity iterations (for each K-inverse solution)
useMG      = 0;
    % 0 ==> CG algorithm preconditioned by MG
    % 1 ==> MG algorithm (no CG); slower but good for comparison
dU_max     = 0.05;
    % Exit iterations if all nodal velocity components change by less than
    % this fraction

% PART 4 - ELEMENT ASSEMBLY
OPTS.method    = 'opt';
    % std = standard element assembly (slow)
    % opt = optimized assembly (blocks of elements; MILAMIN method)
    %       see Dabrowski et al. 2008 (doi:10.1029/2007GC001719)
% when using MILAMIN method:
OPTS.nelblk    = 800;
    % when using 'opt': number of elements assembled at once (block size)
    % Best speed depends on computer architecture (CPU cache)
    % try numbers between 200 and 1000+
OPTS.nel_asmbl = 50000;
    % max number of elements for sparse command
OPTS.nip       = 5;
    % number of integration points
if ~isfield(SETTINGS,'method_eval_visc')
    error(' SETTINGS.method_eval_visc not defined. Check subfunction "calc_element_visc" for choices.');
end
if ~isfield(SETTINGS,'method_eval_dens')
    error(' SETTINGS.method_eval_dens not defined. Check subfunction "calc_element_dens" for choices.');
end
OPTS.method_eval_visc = SETTINGS.method_eval_visc;
OPTS.method_eval_dens = SETTINGS.method_eval_dens;
    % defines the methods used to calculate density and viscosity at each
    % element's integration points
OPTS.return_visc_ip = 0;
    % 1 --> return viscosities at all integration points
OPTS.return_dens_el = 1;
    % 1 --> return mean density in each element

                                                                            tic
% =========================================================================
% ASSEMBLY OF GLOBAL MATRICES
% =========================================================================
[KK,Fb,GG,MM_mu,Fpdil,DensEl,ViscIP] = assembly_stokes_sph_3d...
    (VAR,MESH,PHYSICS,NUMSCALE,OPTS);                                       Uinfo.Tasmble=toc;tic
fprintf(' Assembly of %6i elements in %3i blocks   = %5.1f sec\n',MESH.nel,ceil(MESH.nel/OPTS.nelblk),Uinfo.Tasmble);
KK = KK + tril(KK,-1)'; % Create upper half of symmetric K matrix


% =========================================================================
% FREE SURFACE TERMS
% =========================================================================
switch SETTINGS.top_surface
    case 'free'
        if dt>0
            % Free surface
            alpha    = SETTINGS.fs_alpha; % weighting begin/end of time step
            KKfs_zz  = free_surface_terms(MESH,PHYSICS,NUMSCALE,DensEl,alpha,dt);
            KK       = KK - KKfs_zz; clear KKfs_zz
        end
    case 'fix'
        % Fixed surface

    otherwise
        error(' SETTINGS.top_surface must be either "free" or "fix"');
end


% =========================================================================
% BOUNDARY CONDITIONS
% =========================================================================
nU          = 3*MESH.nnod;
ifix        = UBC.ifix; % list of constrained (fixed) velocity degrees of freedom
Ufix        = UBC.vfix; % list of values at fixed velocity degrees of freedom
ifree       = 1:nU;
ifree(ifix) = []; % list of free velocity degrees of freedom
U           = zeros(nU,1);
U(1:3:end)  = VAR.Ux;   % velocity vector Ux(1), Uy(1), Uz(1), Ux(2), Uy(2), ...
U(2:3:end)  = VAR.Uy;   % velocity vector Ux(1), Uy(1), Uz(1), Ux(2), Uy(2), ...
U(3:3:end)  = VAR.Uz;   % velocity vector Ux(1), Uy(1), Uz(1), Ux(2), Uy(2), ...
U(ifix)     = Ufix;     % prescribed velocity values
P           = VAR.P;    % pressure vector
U0          = U;

% if KK is only lower half:
% KD       = full(diag(KK));
% bc_rhs = KK(:,ifix)*Ufix(:);             % multiply lower part of KK
% bc_rhs = bc_rhs + KK(ifix,:)'*Ufix(:); % multiply "upper" part of KK
% subtract diagonal because it was counted twice
% bc_rhs(ifix) = bc_rhs(ifix) - KD(ifix).*Ufix(:);

% if KK contains both upper and lower triangular part:
bc_rhs = KK(:,ifix)*Ufix(:);


% =========================================================================
% Check if domain is "closed" so that boundary conditions do not allow
% for unconstrained flow in or out of the domain (i.e. zero normal-stress 
% boundary conditions are excluded everywhere).
% If the domain is closed, the pressure is defined up to an arbitrary 
% constant (problem becomes indefinite). During the Conjugate Gradient 
% iterations we will force this constant to be zero by subtracting the mean
% of the pressure solution, which allows CG to solve this 
% positive-indefinite problem.
% =========================================================================
rmPnullspace = is_domain_closed(MESH.PointID,UBC.ifix); %,MESH.GCOORD,MESH.EL2NOD{1});
if rmPnullspace
    rm0space = @(a) a-mean(a);
    fprintf(' BCs DO NOT allow for unconstrained in- or outflow. Will remove nullspace.\n');
else
    fprintf(' BCs DO allow for unconstrained in- or outflow. Will NOT remove nullspace.\n');
end


                                                                            tic
% =========================================================================
% PREPARATION OF STIFFNESS MATRICES ON ALL MULTIGRID LEVELS
% =========================================================================
% Restrict stiffness matrix to all multigrid levels
% Note: KK is converted to a cell array with "nmg" elements; each element
% contains the stiffness matrix on the respective multigrid level.
[KK,Ic2f_U] = restrict_stiffmat(MESH,KK,ifree);                            Uinfo.Trestrict=toc;tic
fprintf(' Restriction of stiffness matrix             = %5.1f sec\n',Uinfo.Trestrict);

% =========================================================================
% CHOLESKY FACTORIZATION OF STIFFNESS MATRIX OON BASE-LEVEL
% =========================================================================
[LbD,permLbD] = cholesky_factorization(KK{end});                            Uinfo.Tchol=toc;tic
fprintf(' Factorization of coarse stiffness matrix    = %5.1f sec\n',Uinfo.Tchol);



% =========================================================================
% PREPARATION OF PRESSURE PRECONDITIONER
% =========================================================================
if PC_Pat==1 || (PC_Pat==4 && OPTS_PC.PC==1)
    DiMM = 1./diag(MM_mu);  % inverse diagonal of MM_mu
end
if PC_Pat==2 || PC_Pat==3 || (PC_Pat==4 && (OPTS_PC.PC==2 || OPTS_PC.PC==3))
    MM_mu = MM_mu + tril(MM_mu,-1)';
    LiMM  = 1./sum(MM_mu,2);  % inverse lumped MM_mu
end



% =========================================================================
%              BEGIN OF VELOCITY/PRESSURE SOLUTION ALGORITHM
%                          (PATERA ALGORITHM)
%
% Important variables used in Preconditioned Conjuate Gradient algorithm
% solving the equation S P = rhs_P
% r     :: residual vector
% d     :: defect vector (approx. to unknown true error of current solution P)
% q     :: search direction in Krylov subspace
% alpha :: step size into search direction q
% beta  :: used to make new search directions q S-orthogonal to all previous q's
% =========================================================================


% =========================================================================
% CALCULATE R.H.S VECTOR
% =========================================================================
Rhs = Fb + GG*P - bc_rhs;
    % buoyancy forces + forces resulting from pressure gradients +
    % velocity boundary conditions
Rhs = Rhs(ifree);


if FigNo
    nsub=2; if PC_Pat==4; nsub=4; end
    figure(FigNo),clf;
end

% =========================================================================
% CALCULATE INITIAL RESIDUAL VECTOR
% =========================================================================
[r,U,time,nit] = pressure_residual(U,Rhs,0); % *NESTED FUNCTION*
nitU    = nit;
time_U  = time;
if rmPnullspace; r = rm0space(r); end
rP_rms  = normdf(r); % norm of residual vector (shown in convergence plot)
tol_Pat = rtol_Pat*rP_rms;
tol_U   = rtol_U*tol_Pat;


% Show profiling data in terminal
% ===============================
% fprintf(' Absolute tolerances to calculate initial rP = %0.2e\n',tol_U0);
% fprintf(' Relative/absolute tolerances for pressure   = %0.2e / %0.2e\n',rtol_Pat,tol_Pat);
% fprintf(' Relative/absolute tolerances for velocity   = %0.2e / %0.2e\n',rtol_U,tol_U);

ind = find(abs(U0)>1e-12*max(abs(U0)));
if isempty(ind)
    dU = max(abs(U));
else
    dU = max(abs(U(ind)./U0(ind)));
end
display_profiling_data(0); clear U0

restart_CG = 1;
                                                                            time_Pat=0;nit_PCPat=0;nit_z=0;profiling=cell(itmax_Pat,1);

% =========================================================================
%  BEGIN OF PRESSURE ITERATIONS
% =========================================================================
for itPat = 1:itmax_Pat
                                                                            tic
    % Pressure preconditioning (approximate error of current P)
    % =========================================================
    [d,z,profiling{itPat}] = precondition_P; % *NESTED FUNCTION*
                                                                            if PC_Pat==4;nit_PCPat=nit_PCPat+profiling{itPat}.nit_PCPat+1;nit_z=nit_z+profiling{itPat}.nit_z;end
    % get new search direction q
    % ==========================
    if restart_CG
        % Initialize Conjugate Gradients
        % ==============================
        q  = d; % define FIRST search direction q
        restart_CG = 0;
        
    else
        rdlast = rd; % save old rd (for calculating beta)
        % Make new search direction q S-orthogonal to all previous q's
        beta = dot(r-rlast,d)/rdlast;  % Polak-Ribiere version 1
%         beta = dot(d-dlast,r)/rdlast;  % Polak-Ribiere version 2
        q    = d + beta*q; % calculate NEW search direction
    end
    
    if rmPnullspace; q = rm0space(q); end
    rd = dot(r,d); % numerator for calculating alpha

    % Perform the S times q multiplication
    % ====================================
    % S cannot be calculated explicitly, since Kinv cannot be formed
    % Sq = S * q = (G' * Kinv * G) * q
    % Hence, the muliplicatipon is done in 3 steps:
    % (1) y   = G*q
    % (2) K z = y (Multigrid-preconditioned Conjugate Gradient solver)
    % (3) Sq  = G'*z
    y  = GG*q;                                                              time_Pat=time_Pat+toc;tic
    [z(ifree),nit,time,rU_rms_vec] = mgpcg_solver...
        (KK,LbD,permLbD,y(ifree),z(ifree),Ic2f_U,MESH.nmg,OPT_mgpcg);       time_U=time_U+time;nitU=nitU+nit;tic
    Sq = GG'*z;

    qSq   = dot(q,Sq); % denominator in calculating alpha
    rlast = r;      % needed for Polak-Ribiere version 1
%     dlast = d;      % needed for Polak-Ribiere version 2
    alpha = rd/qSq; % step size in direction q
    
    % Update solution and residual
    % ============================
    P = P + alpha*q;  % update pressure solution
    U = U + alpha*z;  % update velocity by accumulating alpha*z
    r = r - alpha*Sq; % update residual vector
    if rmPnullspace; r = rm0space(r); end % remove nullspace from both q AND Sq

    rP_rms = normdf(r); % norm of residual vector (shown in convergence plot)
    if FigNo
        figure(FigNo);
        plot(itPat,log10(rP_rms),'k.','Parent',ah1); drawnow;
        plot(log10(rU_rms_vec),'k.-','Parent',ah2); drawnow;
    end
    
    % Show profiling data in terminal
    % ===============================
    display_profiling_data(itPat);
                                                                            time_Pat=time_Pat+toc;
    % Check if solution converged
    % ===========================
    ind = abs(U)>1e-12*max(abs(U));
    dU  = max(abs((alpha*z(ind))./U(ind)));
    if dU<dU_max
        break
    else
        if rP_rms<tol_Pat && itPat>=itmin_Pat
            % =============================================================
            % CALCULATE R.H.S VECTOR
            % =============================================================
            Rhs = Fb + GG*P - bc_rhs;
                % buoyancy forces + forces resulting from pressure gradients
                % + velocity boundary conditions
            Rhs = Rhs(ifree);

            % =============================================================
            % RECALCULATE RESIDUAL VECTOR
            % =============================================================
            tol_Pat  = rtol_Pat*rP_rms;
            tol_U    = rtol_U*tol_Pat;
            [r,U,time,nit] = pressure_residual(U,Rhs,itPat); % *NESTED FUNCTION*
            nitU     = nitU + nit;
            time_U   = time_U + time;
            if rmPnullspace; r = rm0space(r); end
            rP_rms   = normdf(r); % norm of residual vector (shown in convergence plot)
            
            restart_CG  = 1;
            tol_Pat  = rtol_Pat*rP_rms;
            tol_U    = rtol_U*tol_Pat;
        end
    end
end % END OF PRESSURE ITERATIONS

% if rP_rms>tol_Pat
%     fprintf('\nWARNING: Pressure iterations did not converge within %3i iterations: tol=%0.2e rP_rms=%0.2e',...
%             itPat,tol_Pat,rP_rms);
% end

% Check solution
if any(isnan(U)) || any(isnan(P))
	error('Solution containes NaNs !!!!');
end

VAR.Ux = U(1:3:end); % x-velocity solution
VAR.Uy = U(2:3:end); % y-velocity solution
VAR.Uz = U(3:3:end); % z-velocity solution
VAR.P  = P;          % pressure solution

% =========================================================================
% ADJUST MESH FOR SURFACE MOTION
% =========================================================================
if strcmp(SETTINGS.top_surface,'free')
    FigNo      = 88; COMM.nsd = 1;
    [MESH,VAR] = adjust_mesh_fs(VAR,MESH,COMM,dt,SETTINGS,PHYSICS,NUMSCALE,FigNo);
end

fprintf('\n Total time for solving flow problem          = %5.1f sec\n',etime(clock,clock0));
fprintf(' Time spent on pressure problem               = %5.1f sec\n',time_Pat);
fprintf(' Time spent on velocity problem               = %5.1f sec\n',time_U);
fprintf(' Total number of Patera iterations            = %5i\n',itPat);
fprintf(' Total number of velocity iterations          = %5i\n',nitU);
if PC_Pat==4
fprintf(' Total number of inexact Patera iterations    = %5i\n',nit_PCPat-1);
fprintf(' Total number of Z iterations                 = %5i\n',nit_z-1);
time_z=0;time_PCPat=0;
for ii=1:itPat
    time_z     = time_z+profiling{itPat}.time_z;
    time_PCPat = time_PCPat+profiling{itPat}.time_PCPat;
end
fprintf(' Time spent on pressure preconditioning       = %5.1f sec\n',time_z);
fprintf(' Time spent on S-inverse calculation          = %5.1f sec\n',time_PCPat);
end
fprintf(' Total time for solving flow problem          = %6.1f sec\n',etime(clock,clock0));
fprintf('\n');

if FigNo
    saveas(gcf,[SETTINGS.outdir '/Lab01x01_StokesFlowSolverConvergence'],'png');
end

% Save profiling data
Uinfo.tP     = time_Pat;
Uinfo.tU     = time_U;
byte2MB      = 1/1048576;
KK1          = KK{1}; %#ok<NASGU>
tmp          = whos('KK1'); clear KK1
Uinfo.memKK  = tmp.bytes * byte2MB;
tmp          = whos('LbD');
Uinfo.memLL  = tmp.bytes * byte2MB;
tmp          = whos('GG');
Uinfo.memGG  = tmp.bytes * byte2MB;
Uinfo.nel    = MESH.nel;
Uinfo.nPdof  = MESH.nVnod;
Uinfo.nU  = nU;
Uinfo.nitP   = itPat;
Uinfo.rP_rms = rP_rms;
Uinfo.PCPat  = profiling;
Uinfo.nitU   = nitU;

% =========================================================================
%                            NESTED FUNCTIONS
% =========================================================================

function [r,U,time,nit] = pressure_residual(U,Rhs,itPat)
                                                                            
rU_rms_vec=[]; nit=0; time=0;
for it=1:3
    if itPat==0 && it==1
        OPT_mgpcg = struct('useMG',useMG,...
                           'itmax',itmax_U,...
                           'tol'  ,rtol_U,...
                           'is_absolute_tolerance',0);
    else
        OPT_mgpcg = struct('useMG',useMG,...
                           'itmax',itmax_U,...
                           'tol'  ,tol_U,...
                           'is_absolute_tolerance',1);
    end
    
    [U(ifree),nitU0,ttU0,rms_vec] = mgpcg_solver...
        (KK,LbD,permLbD,Rhs,U(ifree),Ic2f_U,MESH.nmg,OPT_mgpcg);           time=time+ttU0;nit=nit+nitU0;rU_rms_vec(end+1:end+length(rms_vec))=rms_vec;tic

    % r = -rhs_P - S p + Fpdil
    %   = -G'*Kinv*Fb - G'*Kinv*G*p + Fpdil , but Fp=G*p and U=Kinv*F (see above)
    % 	= -G'*u_ib - G'*Kinv*Fp + Fpdil     , F=Fb+Fp
    % 	= -G'*U + Fpdil
    r     = -GG'*U + Fpdil; % initial pressure residual vector
    
    if rmPnullspace; r = rm0space(r); end
    rP_rms  = normdf(r); % norm of residual vector (shown in convergence plot)
    tol_Pat = rtol_Pat*rP_rms;
    tol_U   = rtol_U*tol_Pat;
    
    if tol_U>rU_rms_vec(end)
%         fprintf(' Initial pressure residual was calculated.\n');
        break
    end
end

% Plot convergence
% ================
if FigNo
    ah1=subplot(1,nsub,1); plot(itPat,log10(rP_rms),'ro');
    if ~isempty(YLimit); set(gca,'YLim',YLimit); end
    hold all
    set(gca,'FontSize',fs); grid on
    xlabel('Pressure iteration'); ylabel('Norm of pressure residual vector');
    line([0 itmax_Pat],log10([tol_Pat tol_Pat]),...
          'Linestyle','--','Color','r','Linewidth',2);

    ah2=subplot(1,nsub,2); plot(log10(rU_rms_vec),'r.-');
    if ~isempty(YLimit); set(gca,'YLim',YLimit); end
    hold on
    set(gca,'FontSize',fs); grid on
    xlabel('Velocity iteration'); ylabel('Norm of velocity residual vector');
    line([0 itmax_U],log10([tol_U tol_U]),...
          'Linestyle','--','Color','k','Linewidth',2);
end

end % END OF NESTED FUNCTION pressure_residual

% =========================================================================

function [d,z,profiling] = precondition_P
    % Get approximation d to the unknown error of the current
    % pressure solution P (i.e. estimate approx error d from residual r)
    profiling=[]; z = zeros(nU,1);
    switch PC_Pat
        case 0 % No preconditioning (i.e. CG algorithm)
            d = r;
        case 1 % Scaling by diagonal of (1/viscosity)-scaled mass matrix
            d = r.*DiMM;
        case 2 % Scaling by lumped (1/viscosity)-scaled mass matrix
            d = r.*LiMM;
        case 3 % Jacobi iterations on (1/viscosity)-scaled mass matrix
            d  = r .* LiMM;
            for is=1:nit_MM
                Md = MM_mu*d;
                d = d + LiMM .* (r - Md);
            end
        case 4 % inexact Patera algorithm solving S d = r
            OPTS_PC.tol = OPTS_PC.rtol * rP_rms;
            OPTS_PC.tol = max(OPTS_PC.tol,0.75*tol_Pat);
            [d,z,profiling] = inexact_Patera(r,OPTS_PC);
    end
end % END OF NESTED FUNCTION precondition_P

% =========================================================================

function display_profiling_data(flag)
if flag==0
    fprintf(' Profiling of iterative velocity-pressure solution...\n\n');
    if PC_Pat==4
        fprintf('               Patera algorithm                         ||      inexact Patera algorithm     \n');
        fprintf(' -------------------------------------------------------||---------------------------------- \n');
        fprintf('  itP |  time   |    rP    |    rU    |    dU    | nitU ||    rP    | nitP |    rU    | nitU  \n');
        fprintf(' %4i | %7.1f | %0.2e | %0.2e | %0.2e | %4i ||          |      |          |       \n',...
                0,etime(clock,clock0),rP_rms,rU_rms_vec(end),dU,nit);
    else
        fprintf('               Patera algorithm              \n');
        fprintf('  itP |  time   |    rP    |    rU    |    dU    | nitU \n');
        fprintf(' %4i | %7.1f | %0.2e | %0.2e | %0.2e | %4i \n',...
                0,etime(clock,clock0),rP_rms,rU_rms_vec(end),dU,nit);
    end
else
    if PC_Pat==4
        fprintf(' %4i | %7.1f | %0.2e | %0.2e | %0.2e | %4i || %0.2e | %4i | %0.2e | %4i  \n',...
                itPat,etime(clock,clock0),rP_rms,rU_rms_vec(end),dU,nit,...
                profiling{itPat}.rr_rms,profiling{itPat}.nit_PCPat,...
                profiling{itPat}.rZ_rms,profiling{itPat}.nit_z);
        %         1234 | 12345.6 | 1.23e+12 | 1.23e+12 | 1234  || 1.23e+12 | 1234 | 1.23e+12 | 1234
    else
        fprintf(' %4i | %7.1f | %0.2e | %0.2e | %0.2e | %4i \n',...
                itPat,etime(clock,clock0),rP_rms,rU_rms_vec(end),dU,nit);
    end
end
end % END OF NESTED FUNCTION display_profiling_data

% =========================================================================

function [d2,z0,profiling] = inexact_Patera(rr,OPTS_PC)

    % Solves S d = r using Patera's algorithm
    % NOTE:
    %   residual of S d = r  is  rr = r - S * d
    %   with d(:)=0, rr == r
    %   thus input argument r is called rr already
    % NOTE further:
    %   This function is nested inside flow2dD_patera_inexct_mg_p so that
    %   the workspaces are shared. To avoid any conflicts all CG-related
    %   variables that have to remain untouched in the calling funciton
    %   have a suffix "2" in this function (i.e. q is called q2, etc).
    %   These are the critical scalars/vectors:
    %       q      --> q2
    %       beta   --> beta2
    %       rd     --> rd2
    %       rdlast --> rd2last
    %       rlast  --> rlast2
    %       dlast  --> dlast2
    %   This is not required for alpha, Sq, qSq, y or z. These are
    %   recalculated after the call of this function.
    
    nP    = length(r);
    d2    = zeros(nP,1);
    z0    = zeros(nU,1);
    tol_Z = OPTS_PC.rtol_Z*OPTS_PC.tol; % tolerance for inner CG
    
    OPT2_mgpcg                       = OPT_mgpcg;
    OPT2_mgpcg.tol                   = tol_Z;
    OPT2_mgpcg.is_absolute_tolerance = 1;
    
    rr_rms        = normdf(rr);
    rr_rms_vec    = zeros(OPTS_PC.itmax+1,1);
    rr_rms_vec(1) = rr_rms;
    
    if FigNo
        figure(FigNo);
        ah3=subplot(1,nsub,3);
        plot(nit_PCPat,log10(rr_rms),'ro','Markersize',5,'MarkerFaceColor','r'); 
        if itPat==1
            set(gca,'FontSize',fs);
            if ~isempty(YLimit); set(gca,'YLim',YLimit); end
            hold all; grid on
            xlabel('Inexact Patera iteration'); ylabel('Norm of "defect" residual vector');
        end
        set(gca,'XLim',[0 nit_PCPat+OPTS_PC.itmax]);
        line([0 nit_PCPat+OPTS_PC.itmax],log10([OPTS_PC.tol OPTS_PC.tol]),...
             'Linestyle','--','Color','b','Linewidth',2);
        ah4=subplot(1,nsub,4);
        if itPat==1
            if ~isempty(YLimit); set(gca,'YLim',YLimit); end
            hold all; grid on
            set(gca,'FontSize',fs);
            xlabel('Z iteration'); ylabel('Norm of Z residual vector');
        end
        line([0 itmax_U],log10([tol_Z tol_Z]),...
             'Linestyle','--','Color','k','Linewidth',2);
    end
                                                                            time_PCPat=0;itz=0;time_z=0;
    % **********************************
    % Begin of Inexact Patera iterations
    % **********************************
    for it_PCPat = 1:OPTS_PC.itmax
                                                                            tic
        % Pressure preconditioning (approximate error of current P)
        % =========================================================
        dd = precondition_d; % *NESTED FUNCTION*

        % get new search direction q2
        % ===========================
        if it_PCPat==1
            % Initialize Conjugate Gradients
            % ==============================
            q2 = dd; % define FIRST search direction q

        else
            rd2last = rd2; % save old rd (for calculating beta)
            % Make new search direction q S-orthogonal to all previous q's
            beta2 = dot(rr-rlast2,dd)/rd2last;  % Polak-Ribiere version 1
%             beta2 = dot(dd-dlast2,rr)/rd2last; % Polak-Ribiere version 2
            q2    = dd + beta2*q2; % calculate NEW search direction
        end

        if rmPnullspace; q2 = rm0space(q2); end
        rd2 = dot(rr,dd); % numerator for calculating alpha

        % Perform the S times q multiplication
        % ====================================
        % S cannot be calculated explicitly, since Kinv cannot be formed
        % Sq = S * q = (G' * Kinv * G) * q
        % Hence, the muliplicatipon is done in 3 steps:
        % (1) y   = G*q
        % (2) K z = y (Multigrid-preconditioned Conjugate Gradient solver)
        % (3) Sq  = G'*z
        y  = GG*q2;                                                         time_PCPat=time_PCPat+toc;tic
        z  = zeros(nU,1);
        [z(ifree),itz1,ttz1,rz_rms_vec] = mgpcg_solver...
            (KK,LbD,permLbD,y(ifree),z(ifree),Ic2f_U,MESH.nmg,OPT2_mgpcg);time_z=time_z+ttz1;itz=itz+itz1;tic
        Sq = GG'*z;

        qSq    = dot(q2,Sq); % denominator in calculating alpha
        rlast2 = rr;      % needed for Polak-Ribiere version 1
%         dlast2 = d2;      % needed for Polak-Ribiere version 2
        alpha  = rd2/qSq; % step size in direction q

        % Update solution and residual
        % ============================
        d2 = d2 + alpha*q2; % update pressure "defect"
        rr = rr - alpha*Sq; % update residual of "S d = r"
        z0 = z0 + alpha*z;  % update velocity correction by accumulating alpha*z
        if rmPnullspace; rr = rm0space(rr); end % remove nullspace from both q AND Sq

        rr_rms = normdf(rr); % norm of residual vector (shown in convergence plot)
        rr_rms_vec(it_PCPat+1) = rr_rms;        
                                                                            time_PCPat=time_PCPat+toc;
        if FigNo
            plot(nit_PCPat+it_PCPat,log10(rr_rms),'b.','Parent',ah3); drawnow;
            plot(log10(rz_rms_vec),'k.-','Parent',ah4); drawnow;
        end
        % Check if solution converged
        % ===========================
        if rr_rms<OPTS_PC.tol && it_PCPat>=OPTS_PC.itmin
            break
        end
    end % End of Inexact Patera iterations
    % ************************************

    if rr_rms>OPTS_PC.tol
        fprintf('WARNING: Inexact Patera iterations did not converge within %3i iterations: tol=%0.2e rr_rms=%0.2e\n',...
                it_PCPat,OPTS_PC.tol,rr_rms);
    else
%         fprintf(' Total time for solving flow problem = %5.1f sec\n',etime(clock,clock0)); %...
%         fprintf(' Total number of Patera iterations   = %5i\n',it_PCPat);
%         fprintf(' Total number of velocity iterations = %5i\n',itz);
    end
    
    profiling.tol_PCPat  = OPTS_PC.tol;
    profiling.tol_Z      = tol_Z;
    profiling.rZ_rms     = rz_rms_vec(end);
    profiling.rr_rms     = rr_rms;
    profiling.rr_rms_vec = rr_rms_vec(1:it_PCPat+1);
    profiling.nit_PCPat  = it_PCPat;
    profiling.time_PCPat = time_PCPat;
    profiling.nit_z      = itz;
    profiling.time_z     = time_z;
    
    %%%%%%%%%%%%%%%% NESTED FUNCTIONS %%%%%%%%%%%%%%%%%

    function dd = precondition_d
        % Get approximation d to the unknown error of the current
        % pressure solution P (i.e. estimate approx error d from residual r)
        switch OPTS_PC.PC
            case 0 % No preconditioning (i.e. CG algorithm)
                dd = rr;
            case 1 % Scaling by diagonal of (1/viscosity)-scaled mass matrix
                dd = rr.*DiMM;
            case 2 % Scaling by lumped (1/viscosity)-scaled mass matrix
                dd = rr.*LiMM;
            case 3 % Jacobi iterations on (1/viscosity)-scaled mass matrix
                dd = rr .* LiMM;
                for is=1:nit_MM
                    Md = MM_mu*dd;
                    dd = dd + LiMM .* (rr - Md);
                end
        end
    end % END OF NESTED FUNCTION precondition_d
    
end % END OF NESTED FUNCTION inexact_Patera

end % END OF FUNCTION flow3d_patera_mg_inexact