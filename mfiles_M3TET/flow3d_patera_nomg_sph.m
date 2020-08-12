function [VAR,Uinfo] = flow3d_patera_nomg_sph(VAR,MESH,UBC,SETTINGS,PHYSICS,NUMSCALE,dt)
% Usage: [VAR,Uinfo] = flow3d_patera_nomg_sph(VAR,MESH,UBC,SETTINGS,PHYSICS,NUMSCALE,dt)
%
% Purpose: Stokes flow solver (calculate velocity and pressure solution)
%
% Input:
%   VAR      : [structure] : structure containing all major variables
%   MESH     : [structure] : FE mesh parameters
%   UBC      : [structure] : velocity boundary conditions
%   SETTINGS : [structure] : model parameters
%   PHYSICS  : [structure] : physical properties
%   NUMSCALE : [structure] : numerical scaling parameters
%   dt       : [double]    : time step (not used yet)
%
% Output:
%   VAR      : [structure] : structure containing all major variables
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
% - velocity solution: Conjugate Gradient algorithm, direct solver, or 
%                      algebratic multigrid (AGMG)
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
% JH/JPM Feb 2011
% JH Mar 2011
% JH May 2014: - added Yvan Notay's algebraic multigrid solver AGMG as
%                alternative solver
% JH Feb 2016: - fine tuning of tolerances
%              - verified convergence for zero velocities at CMB and surface
%
                                                                            clock0=clock;tic
% =========================================================================
% SETUP FOR STOKES FLOW SOLVER
% =========================================================================
FigNo    = 144;
    % number of figure showing convergence; ==0 --> no plot
fs       = 10;
    % font size in convergence plot
YLimit   = [-4 10];
    % ylimits for convergence plots
    
% PART 1 - PRESSURE ITERATIONS
rtol_Pat   = 1e-4;
    % RELATIVE tolerance for pressure solution
    % (relative to norm of initial pressure residual)
pc_Pat     = 3;
    % 0 ==> no preconditioner (d=r)
    % 1 ==> d = r./diag(MM);  scaling by diagonal of mass matrix
    % 2 ==> d = r./LMM; scaling by lumped mass matrix; LMM=sum(MM,2)
    % 3 ==> d = d + (r - MM*d)./LMM; Jacobi iterations on MM 
    % 4 ==> zero-infill Cholesky factorization of MM
    % 5 ==> incomplete Cholesky factorization of MM
    % 6 ==> Cholesky factorization of MM
% wght_MM   = [1];
% wght_MM   = [0.25 0.4 2/3];
wght_MM   = [0.25 0.5 0.75 1 0.75 0.5 0.25];
    % number of mass matrix iterations
if pc_Pat==4
    itmin_Pat = 1;   % minimum number of pressure iterations
    itmax_Pat = 20;  % maximum number of pressure iterations
else
    itmin_Pat = 3;   % minimum number of pressure iterations
    itmax_Pat = 120; % maximum number of pressure iterations
end

% PART 2 - VELOCITY SOLUTION
solver_U = 'chol';
    % 'bslash' = Matlab's backslash
    % 'chol' = Cholesky direct solver (only for VERY SMALL MESHES!!!)
    % 'pcg'  = Preconditoned Conjugate Gradient solver (you must specify a 
    %          preconditioner "pc_Vel")
    % 'agmg' = algebraic multigrid (NOT YET WORKING)
pc_Vel     = '0chol'; % only used if solver_U==pcg
    % '0chol': zero-infill incomplete Cholesky factorization
    % 'ichol': drop level incomplete Cholesky factorization
    % 'symGS': symmetric Gauss-Seidel
    % 'diag' : diagonal scaling (Jacobi preconditioner)
    % 'noPC' : no preconditioning (CG algorithm)
rtol_U   = 5e-3;
    % tolerance for velocity solution (RELATIVE to pressure tolerance; use ~1e-2)
itmax_U  = 100;
    % maximum number of iterations for pcg and agmg
    
% PART 3 - ELEMENT ASSEMBLY
OPTS.method    = 'opt';
    % std = standard element assembly (slow)
    % opt = optimized assembly (blocks of elements; MILAMIN method)
    %       see Dabrowski et al. 2008 (doi:10.1029/2007GC001719)
% when using MILAMIN method:
OPTS.nelblk    = 1400;
    % when using 'opt': number of elements assembled at once (block size)
    % Best speed depends on computer architecture (CPU cache)
    % try numbers between 200 and 1000+
OPTS.use_mutils = 1;
    % 1 --> use MUTILS' "sparse_create"
COMM          = init_COMM();
OPTS.nthreads = COMM.nthreads; clear COMM
    % number of threads used by mutils
OPTS.nip        = 14;
    % number of integration points (4, 5, 8, 10, 11, 14)
OPTS.type_Dmat  = '23rd'; % '221' or '23rd'
    % D-matrix formulation (extraction of deviatoric strain rates)
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
OPTS.return_dens_el = 0;
    % 1 --> return mean density in each element
fidl = SETTINGS.fid_log;
if strcmp(SETTINGS.show_Figs,'no')
    FigNo = 0;
end

% Set to 'yes' to make gravity point into same direction within each
% element (only works with OPTS.nip==5 !!!!)
OPTS.const_Fg_per_el = 'no';

                                                                            tic
%==========================================================================
% ASSEMBLY OF GLOBAL MATRICES
%==========================================================================
switch SETTINGS.jacobian
    case 'standard'
        [KK,Fb,GG,MM_mu,Fpdil] = assembly_stokes_sph_3d...
            (VAR,MESH,PHYSICS,NUMSCALE,OPTS);                              Uinfo.Tasmble=toc;tic
        fprintf(fidl,' Assembly of %6i elements in %3i blocks   = %5.1f sec\n',...
            MESH.nel,ceil(MESH.nel/OPTS.nelblk),Uinfo.Tasmble);
    case 'double'
        [KK,Fb,GG,MM_mu,Fpdil] = assembly_stokes_sph_3d_double_jacobian...
            (VAR,MESH,PHYSICS,NUMSCALE,OPTS);                              Uinfo.Tasmble=toc;tic
        fprintf(fidl,' Assembly of %6i elements in %3i blocks   = %5.1f sec\n',...
            MESH.nel,ceil(MESH.nel/OPTS.nelblk),Uinfo.Tasmble);
end

if strcmp(SETTINGS.plate_model,'WJM')
    %======================================================================
    % PROBLEM SIZE, ALLOCATE SOLUTION VECTORS
    %======================================================================
    nU             = 3*MESH.nnod;
    U_xyz          = zeros(nU,1);
    U_xyz(1:3:end) = VAR.Ux; % velocity vector Ux(1), Uy(1), Uz(1), Ux(2), Uy(2), ...
    U_xyz(2:3:end) = VAR.Uy; % velocity vector Ux(1), Uy(1), Uz(1), Ux(2), Uy(2), ...
    U_xyz(3:3:end) = VAR.Uz; % velocity vector Ux(1), Uy(1), Uz(1), Ux(2), Uy(2), ...
    P              = VAR.P;  % pressure vector
    %======================================================================
    % MODIFY MATRICES AND FORCE VECTOR USING THE ROTATION MATRIX
    %======================================================================
    RR      = MESH.RR_xyz2rEN;
    % RR      = speye(size(RR));
    is_tril = full(all(KK(1,2:end)==0));
    KK      = KK + tril(KK,-1)'; % Create upper half of symmetric K matrix
    KK      = RR*KK*RR';
    if is_tril
        KK = tril(KK);
    end
    GG = RR * GG;
    Fb = RR * Fb;
    U  = RR * U_xyz;
elseif sum(strcmp(SETTINGS.plate_model,{'GPlates','Debug'}))
    %======================================================================
    % PROBLEM SIZE, ALLOCATE SOLUTION VECTORS
    %======================================================================
    nU             = 3*MESH.nnod;
    U              = zeros(nU,1);
    P              = VAR.P;  % pressure vector
    % MODIFY MATRICES AND FORCE VECTOR USING THE ROTATION MATRIX
    RR      = MESH.RR;
    % RR      = speye(size(RR));
    is_tril = full(all(KK(1,2:end)==0));
    KK      = KK + tril(KK,-1)'; % Create upper half of symmetric K matrix
    KK      = RR'*KK*RR;
    if is_tril
        KK = tril(KK);
    end
    GG = RR' * GG;
    Fb = RR' * Fb;
end
                                                                            Uinfo.Trotate=toc;
fprintf(fidl,' Rotation of stiffness matrix                = %5.1f sec\n',Uinfo.Trotate);

%==========================================================================
% Check if domain is "closed" so that boundary conditions do not allow
% for unconstrained flow in or out of the domain (i.e. zero normal-stress 
% boundary conditions are excluded everywhere).
% If the domain is closed, the pressure is defined up to an arbitrary 
% constant (problem becomes indefinite). During the Conjugate Gradient 
% iterations we will force this constant to be zero by subtracting the mean
% of the pressure solution, which allows CG to solve this 
% positive-indefinite problem.
%==========================================================================
rmPnullspace = 1; % IN A SPHERE THE DOMAIN IS CLOSED
if rmPnullspace
    rm0space = @(a) a-mean(a);
    fprintf(fidl,' BCs DO NOT allow for unconstrained in- or outflow. Will remove nullspace.\n');
else
    fprintf(fidl,' BCs DO allow for unconstrained in- or outflow. Will NOT remove nullspace.\n');
end

%==========================================================================
% BOUNDARY CONDITIONS
%==========================================================================
ifix        = UBC.ifix; % list of constrained (fixed) velocity degrees of freedom
Ufix        = UBC.vfix; % list of values at fixed velocity degrees of freedom
ifree       = 1:nU;
ifree(ifix) = [];       % list of free velocity degrees of freedom
U(ifix)     = Ufix;     % prescribed velocity values

% Subtract boundary condition terms from r.h.s.:
% (if K1 would contain both upper and lower triangular parts it would 
% simply be: bc_rhs = K1(:,ifix)*Ufix(:);)
KD           = full(diag(KK));
bc_rhs       = KK(:,ifix)*Ufix(:);           % multiply lower part of KK
bc_rhs       = bc_rhs + KK(ifix,:)'*Ufix(:); % multiply "upper" part of KK
bc_rhs(ifix) = bc_rhs(ifix) - KD(ifix).*Ufix(:);% subtract diagonal (counted twice)

%==========================================================================
% CALCULATE R.H.S VECTOR
%==========================================================================
Rhs = Fb + GG*P - bc_rhs;
    % buoyancy forces + forces resulting from pressure gradients +
    % velocity boundary conditions
    
Rhs = Rhs(ifree);
KK  = KK(ifree,ifree);

if strcmp(solver_U,'pcg')
    opts_spc.symmetric = 1;
    opts_spc.nthreads  = 1;
    KK = sparse_convert(KK,opts_spc);
elseif ~strcmp(solver_U,'chol')
    KK = KK + tril(KK,-1)'; % Create upper half of symmetric K matrix
end


%==========================================================================
% PREPARATION OF PRESSURE PRECONDITIONER
%==========================================================================
switch pc_Pat
    case 1
        DiMM  = 1./diag(MM_mu); % inverse diagonal of MM_mu
    case 2
        MM_mu = MM_mu + tril(MM_mu,-1)'; % create upper half
        LiMM  = 1./sum(MM_mu,2); % inverse lumped MM
    case 3
        DiMM  = 1./diag(MM_mu); % inverse diagonal of MM_mu
        MM_mu = MM_mu + tril(MM_mu,-1)'; % create upper half
    case 4
        % zero-infill
        OPTS_ichol.type = 'nofill';
        perm_MM = amd(MM_mu);
        LL_mu   = ichol(MM_mu(perm_MM,perm_MM),OPTS_ichol);
    case 5
        % incomplete with droplevel tolerance
        OPTS_ichol.type    = 'ict';
        OPTS_ichol.droptol = 1e-8*max(abs(diag(MM_mu)));
        perm_MM = amd(MM_mu);
        LL_mu   = ichol(MM_mu(perm_MM,perm_MM),OPTS_ichol);
    case 6
        [LL_mu,perm_MM] = cholesky_factorization(MM_mu);
end


%==========================================================================
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
%==========================================================================

if FigNo
    figure(FigNo); clf;
end

%==========================================================================
% CALCULATE INITIAL RESIDUAL VECTOR
%==========================================================================
[r,U,time,nit] = pressure_residual(U,Rhs,0); % *NESTED FUNCTION*
nitU    = nit;
time_U  = time;
if rmPnullspace; r = rm0space(r); end
norm_rP = normdf(r); % norm of residual vector (shown in convergence plot)
tol_Pat = rtol_Pat*norm_rP;
tol_U   = rtol_U*tol_Pat;


% Show profiling data in terminal
% ===============================
% fprintf(' Final tol. for calculating pressure residual = %0.2e\n',rtol_U);
% fprintf(' Rel./abs. tolerances for pressure            = %0.2e / %0.2e\n',rtol_Pat,tol_Pat);
% fprintf(' Rel./abs. tolerances for velocity            = %0.2e / %0.2e\n',rtol_U*rtol_Pat,tol_U);


dU = max(abs(U));
display_profiling_data(0);

restart_CG = 1;
                                                                            time_Pat=0;
%==========================================================================
%  BEGIN OF PRESSURE ITERATIONS
%==========================================================================
for itPat = 1:itmax_Pat
                                                                            tic
    % Pressure preconditioning (approximate error of current P)
    % =========================================================
    d  = precondition_P(r); % *NESTED FUNCTION*
    
    % get new search direction q
    % ==========================
    if restart_CG
        % Initialize Conjugate Gradients
        % ==============================
        q = d; % define FIRST search direction q
        restart_CG = 0;
        
    else
        rdlast = rd; % save old rd (for calculating beta)
        % Make new search direction q S-orthogonal to all previous q's
        beta   = dot(r-rlast,d)/rdlast;  % Polak-Ribiere version 1
    %     beta   = dot(d-dlast,r)/rdlast;  % Polak-Ribiere version 2
        q      = d + beta*q; % calculate NEW search direction
    end
    
    if rmPnullspace; q = rm0space(q); end
    rd = dot(r,d); % numerator for calculating alpha

    % Perform the S times q multiplication
    % ====================================
    % S cannot be calculated explicitly, since Kinv cannot be formed
    % Sq = S * q = (G' * Kinv * G) * q
    % Hence, the muliplicatipon is done in 3 steps:
    % (1) y   = G*q
    % (2) K z = y (Cholesky direct solver)
    % (3) Sq  = G'*z
    y         = GG*q;                                                       time_Pat=time_Pat+toc;tic
    z         = zeros(nU,1);
    switch solver_U
        case 'chol'
            z(ifree) = cholesky_solve(y(ifree),LL,perm);             time_U=toc;
        case 'pcg'
            OPT_pcg.tol        = tol_U; % tolerance for CG solution
            OPT_pcg.is_abs_tol = 1;     % absolute tolerance
            [z(ifree),profiling] = pcg_solver...
                (KK,y(ifree),z(ifree),OPT_pcg);time_U=profiling.time(end);  nitU=nitU+profiling.nit;
            rU_rms_vec = profiling.rrms;
            nit        = profiling.nit;
        case 'agmg'
            icg     = 1;
            rtol_U  = 1e-8;
            verbose = 0;
            [z(ifree),flag,relres,nit,rU_rms_vec] = agmg...
                (KK,y(ifree),icg,rtol_U,itmax_U,verbose);                   time_U=time_U+toc;nitU=nitU+nit;
        case 'bslash'
            z(ifree) = KK \ y(ifree);                                       time_U=toc;
    end

    Sq    = GG'*z;
    qSq   = dot(q,Sq); % denominator in calculating alpha
    rlast = r;         % needed for Polak-Ribiere version 1
%     dlast = d;         % needed for Polak-Ribiere version 2
    alpha = rd/qSq;    % steps size in direction q
    
    % Update solution and residual
    % ============================
    P = P + alpha*q;  % update pressure solution
    U = U + alpha*z;  % update velocity by accumulating alpha*z
    r = r - alpha*Sq; % update residual vector
    if rmPnullspace; r = rm0space(r); end % remove nullspace from both q AND Sq
                                                                            time_Pat=time_Pat+toc;
    norm_rP = normdf(r); % norm of residual vector (shown in convergence plot)
    if FigNo
        figure(FigNo);
        plot(itPat,log10(norm_rP),'k.','Parent',ah1); drawnow;
        if ~ismember(solver_U,{'chol','bslash'})
            semilogy(log10(rU_rms_vec),'k.-','Parent',ah2); drawnow;
        end
    end
    
    % Show profiling data in terminal
    % ===============================
    dU = max(abs(alpha*z));
    display_profiling_data(itPat);
    
    % Check if solution converged
    % ===========================
    if norm_rP<tol_Pat && itPat>=itmin_Pat
       break
    end
end % END OF PRESSURE ITERATIONS

% if norm_rP>tol_Pat
%     fprintf('\nWARNING: Pressure iterations did not converge within %3i iterations: tol=%0.2e norm_rP=%0.2e',...
%             itPat,tol_Pat,norm_rP);
% end

% Check solution
if any(isnan(U)) || any(isnan(P))
	error('Solution containes NaNs !!!!');
end
if strcmp(SETTINGS.plate_model,'WJM')
    VAR.Ur = U(1:3:end); % Radial component
    VAR.Ue = U(2:3:end); % East component
    VAR.Un = U(3:3:end); % North component
    U_xyz  = RR'*U;      % rotate velocity solution back to xyz-coordinates
elseif sum(strcmp(SETTINGS.plate_model,{'GPlates','Debug'}))
    VAR.Uth = U(1:3:end); % theta component
    VAR.Uph = U(2:3:end); % phi component
    VAR.Ur  = U(3:3:end); % radial component
    U_xyz   = RR*U;       % rotate velocity solution back to xyz-coordinates
end
VAR.Ux = U_xyz(1:3:end); % x-velocity solution
VAR.Uy = U_xyz(2:3:end); % y-velocity solution
VAR.Uz = U_xyz(3:3:end); % z-velocity solution
VAR.P  = P;          % pressure solution

% % Add lithostatic pressure
% VAR.P  = VAR.P + (PHYSICS.g * MESH.GCOORD(3,1:nP)' * NUMSCALE.Dens0);
% % Shift pressure so that P==0 at the top
% VAR.P  = VAR.P - min(VAR.P);

fprintf(fidl,'\n Total time for solving flow problem          = %5.1f sec\n',etime(clock,clock0));
fprintf(fidl,' Time spent on pressure problem               = %5.1f sec\n',time_Pat);
fprintf(fidl,' Time spent on velocity problem               = %5.1f sec\n',time_U);
fprintf(fidl,' Total number of Patera iterations            = %5i\n',itPat);
fprintf(fidl,' Total number of velocity iterations          = %5i\n',nitU);
if pc_Pat==4
fprintf(fidl,' Total number of inexact Patera iterations    = %5i\n',nit_PCPat-1);
fprintf(fidl,' Total number of Z iterations                 = %5i\n',nit_z-1);
time_z=0;time_PCPat=0;
for ii=1:itPat
    time_z     = time_z+profiling{itPat}.time_z;
    time_PCPat = time_PCPat+profiling{itPat}.time_PCPat;
end
fprintf(fidl,' Time spent on pressure preconditioning       = %5.1f sec\n',time_z);
fprintf(fidl,' Time spent on S-inverse calculation          = %5.1f sec\n',time_PCPat);
end
fprintf(fidl,'\n');

if FigNo
    saveas(gcf,[SETTINGS.outdir '/Lab01x01_StokesFlowSolverConvergence'],'png');
end



% Save profiling data
% *******************
Uinfo.tP     = time_Pat;
Uinfo.tU     = time_U;
byte2MB      = 1/1048576;
tmp          = whos('KK');
Uinfo.memKK  = tmp.bytes * byte2MB;
if strcmp(solver_U,'chol')
    tmp         = whos('LL');
    Uinfo.memLL = tmp.bytes * byte2MB;
end
Uinfo.memGG  = tmp.bytes * byte2MB;
Uinfo.nel    = MESH.nel;
Uinfo.nPdof  = MESH.nVnod;
Uinfo.nUdof  = nU;
Uinfo.nitP   = itPat;
Uinfo.norm_rP = norm_rP;

%==========================================================================
%                            NESTED FUNCTIONS
%==========================================================================

function [r,U,time,nit] = pressure_residual(U,Rhs,itPat)

rU_rms_vec=[]; nit=0; time = 0;
for it=1:3     
    if strcmp(solver_U,'chol')
        if nU>650000
            error(' The Cholesky direct solver must only be used for small 3D problems.');
        end
                                                                                tic
        [LL,perm] = cholesky_factorization(KK);
        fprintf(fidl,' Factorization of stiffness matrix            = %5.1f sec\n',toc);
        U(ifree) = cholesky_solve(Rhs,LL,perm);                             time=time+toc;nit=nit+1;
        rU_rms_vec = 0;
    elseif strcmp(solver_U,'bslash')
                                                                            tic
        U(ifree) = KK \ Rhs;                                                time=time+toc;nit=nit+1;
        rU_rms_vec = 0;
    else
        switch solver_U
            case 'pcg'
                if itPat==0 && it==1
                    OPT_pcg.is_abs_tol = 0;
                    OPT_pcg.tol = rtol_U; % tolerance for CG solution
                else
                    OPT_pcg.is_abs_tol = 1;
                    OPT_pcg.tol = tol_U; % tolerance for CG solution
                end
                OPT_pcg.itmax = itmax_U; % maximum number of CG iterations
                OPT_pcg.PC    = pc_Vel; % preconditioner 
                [U(ifree),profiling] = pcg_solver(KK,Rhs,U(ifree),OPT_pcg);     time=time+profiling.time(end);nit=nit+profiling.nit;
                rms_vec  = profiling.rrms;
                
            case 'agmg'
                icg     = 10;                                                    tic
                    % 1, use CG and performs some simplifications based on
                    % the assumption that the coefficient matrix A is symmetric;
                    % should be used only when A is symmetric and positive definite.
                    % >=2, use GCR restarted each RESTART iterations.
                    % 0 or [], then agmg use the default, which is
                    % GCR restarted each 10 iterations.
                verbose = 0;
                if itPat==0 && it==1
                    tol = 0.01*rtol_U; % tolerance for CG solution
                else
                    tol = 0.9*rtol_U*tol_U/rms_vec(end); % tolerance for CG solution
                end
                [U(ifree),flag,relres,nit0,rms_vec] = agmg...
                    (KK,Rhs,icg,tol,itmax_U,verbose,U(ifree));              time=time+toc;nit=nit+nit0;
        end
                                                                            rU_rms_vec(end+1:end+length(rms_vec))=rms_vec;
    end

    % r = -rhs_P - S p + Fpdil
    %   = -G'*Kinv*Fb - G'*Kinv*G*p + Fpdil , but Fp=G*p and U=Kinv*F (see above)
    % 	= -G'*u_ib - G'*Kinv*Fp + Fpdil     , F=Fb+Fp
    % 	= -G'*U + Fpdil
    r     = -GG'*U + Fpdil; % initial pressure residual vector

    if rmPnullspace; r = rm0space(r); end
    norm_rP = normdf(r); % norm of residual vector (shown in convergence plot)
    tol_Pat = rtol_Pat*norm_rP;
    tol_U   = rtol_U*tol_Pat;

    if tol_U>rU_rms_vec(end)
%         fprintf(' Initial pressure residual was calculated.\n');
        break
    end
end

% Plot convergence
% ================
if FigNo
    if ~ismember(solver_U,{'chol','bslash'})
        ah1=subplot(1,2,1);
    end
    plot(itPat,log10(norm_rP),'ro');
    hold all
    set(gca,'FontSize',fs); grid on
    xlabel('Pressure iteration'); ylabel('Norm of pressure residual vector');
    line([0 itmax_Pat],log10([tol_Pat tol_Pat]),...
          'Linestyle','--','Color','r','Linewidth',2);
    if ~isempty(YLimit);
        set(gca,'YLim',YLimit);
        if YLimit(1) > log10(tol_Pat)
            YLimit(1) = floor(log10(tol_Pat));
            set(gca,'YLim',YLimit);
        end
    end

    if ismember(solver_U,{'chol','bslash'})
        ah1 = gca;
    else
        ah2=subplot(1,2,2); plot(log10(rU_rms_vec),'r.-');
        hold on
        set(gca,'Xlim',[0 itmax_Pat],'FontSize',fs); grid on
        if ~isempty(YLimit); set(gca,'YLim',YLimit); end
        xlabel('Velocity iteration'); ylabel('Norm of velocity residual vector');
        line([0 itmax_U],log10([tol_U tol_U]),...
              'Linestyle','--','Color','k','Linewidth',2);
    end
end

end % END OF NESTED FUNCTION pressure_residual

% =========================================================================

function d = precondition_P(r)
    % Get approximation d to the unknown error of the current
    % pressure solution P (i.e. estimate approx error d from residual r)
    switch pc_Pat
        case 0 % no preconditioning (i.e. CG algorithm)
            d = r;
        case 1 % scaling by diagonal of (1/viscosity)-scaled mass matrix
            d = r.*DiMM;
        case 2 % scaling by lumped (1/viscosity)-scaled mass matrix
            d = r.*LiMM;
        case 3 % Jacobi iterations on (1/viscosity)-scaled mass matrix
            d  = r .* DiMM;
            for is=1:length(wght_MM)
                d  = d + wght_MM(is) * DiMM .* (r - MM_mu*d);
            end
        case 7
            d  = r ./ diag(SS);
            ns = 5;
            for is=1:ns
                d = d + (r - SS*d)./diag(SS);
            end
%             d = cholesky_substitution(r,SS,perm_S);
        otherwise
            d = cholesky_substitution(r,LL_mu,perm_MM);
    end
end % END OF NESTED FUNCTION precondition_P

% =========================================================================

function display_profiling_data(flag)
    
if strcmp(solver_U,'chol') || strcmp(solver_U,'bslash')
    rU_rms = 0;
else
    rU_rms = rU_rms_vec(end);
end

if flag==0
    fprintf(fidl,' Profiling of iterative velocity-pressure solution...\n\n');
    fprintf(fidl,'               Patera algorithm              \n');
    fprintf(fidl,'  itP |  time   |    rP    |    rU    |    dU    | nitU \n');
    fprintf(fidl,' %4i | %7.1f | %0.2e | %0.2e | %0.2e | %4i \n',...
            0,etime(clock,clock0),norm_rP,rU_rms,dU,nit);
else
    fprintf(fidl,' %4i | %7.1f | %0.2e | %0.2e | %0.2e | %4i \n',...
            itPat,etime(clock,clock0),norm_rP,rU_rms,dU,nit);
end
end % END OF NESTED FUNCTION display_profiling_data

end % END OF FUNCTION flow3d_patera_nomg