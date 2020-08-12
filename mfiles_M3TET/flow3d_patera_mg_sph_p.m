function [VAR,Uinfo] = flow3d_patera_mg_sph_p(VAR,MESH,COMM,UBC,SETTINGS,PHYSICS,NUMSCALE,dt)
% Usage: [VAR,Uinfo] = flow3d_patera_mg_sph_p(VAR,MESH,COMM,UBC,SETTINGS,PHYSICS,NUMSCALE,dt)
%
% Purpose: Stokes flow solver (calculate velocity and pressure solution)
%
% Input:
%   VAR      : [structure] : structure containing all major variables
%   MESH     : [structure] : FE mesh parameters
%   COMM     : [structure] : inter-subdomain communication data
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
% - parallel version
% - elements: 3D Taylor-Hood triangles (10 velocity nodes, 4 pressure nodes)
% - velocity solution: Conjugate Gradient algorithm, preconditioned by a
%                      single V-cylce of geometric multigrid
% - pressure solution: Preconditioned Conjugate Gradient (PCG) solver
%   (for available preconditioners see "setup"-block below)
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
% JH Mar/Apr 2011
% JH May 2014: - max. velocity change during iterations now criterium for
%                ending iterations
%              - added free surface (algorithm by Miguel & Jason)
% JH May 2014: - included sparse-matrix-vector-multiplication from MUTILS
%              - added Yvan Notay's algebraic multigrid solver AGMG as
%                alternative solver on the coarsest mesh (it is slower than
%                Cholesky if the mesh is too small)
% JH Feb 2016: - fine tuning of tolerances
%              - verified convergence for zero velocities at CMB and surface
%
                                                                            clock0=clock;tic
% =========================================================================
% SETUP FOR STOKES FLOW SOLVER
% =========================================================================
FigNo    = 56;
    % number of figure showing convergence; ==0 --> no plot
fs       = 10;
    % font size in convergence plot
YLimit   = [-4 10];
    % y-axis limits in convergence plot

% PART 1 - PRESSURE ITERATIONS
rtol_Pat   = 1e-4;
    % RELATIVE tolerance for pressure solution
    % (relative to norm of initial pressure residual)
% rtol_restart = 5e-2;
rtol_restart = 0.8*rtol_Pat;
    % tolerance after which outer CG is restarted
pc_Pat     = 3;
    % 0 ==> no preconditioner (d=r)
    % 1 ==> d = r./diag(MM);  scaling by diagonal of mass matrix
    % 2 ==> d = r./LMM; scaling by lumped mass matrix; LMM=sum(MM,2)
    % 3 ==> d = d + (r - MM*d)./LMM; Jacobi iterations on MM
wght_MM   = [1];
% wght_MM   = [2/3 1 2/3];
% wght_MM   = [0.25 0.4 2/3];
% wght_MM   = [0.25 0.5 0.75 1 0.75 0.5 0.25];
    % number of mass matrix iterations when using pc_Pat==3 or pc_Pat==4
itmin_Pat = 3;   % minimum number of pressure iterations
itmax_Pat = 120; % maximum number of pressure iterations

% PART 2 - VELOCITY SOLUTION
rtol_U     = 5e-3;
    % tolerance for velocity solution (RELATIVE to pressure tolerance;
    % use ~0.1; will be adjusted if not low enough)
itmax_U    = 20;
    % maximum number of velocity iterations (for each K-inverse solution)
useMG      = 0;
    % 0 ==> CG algorithm preconditioned by MG
    % 1 ==> MG algorithm (no CG); slower but good for comparison
coarse_solver = 'chol'; % use "chol" !!!!!!!!!
    % Choose solver on coarsest multigrid level
    % 'chol'  --> Cholesky direct solver
    % 'ichol' --> incomplete Cholesky factorization
    % 'agmg'  --> Yvan Notay's algebraic multigrid solver AGMG
    % 'bslash' --> Matlab's backslash (use only for testing)
    % 'jacobi' --> Jacobi iterations
use_mutils = 1;
    % 1 --> use MUTILS tools for better preformance
    % 0 --> use Matlab's functions (slower)

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
OPTS.nel_asmbl = 1e8;
    % max number of elements for sparse command
    % use ~50000 when using sparse or sparse2 assembly
    % use 1e8 when using MUTILS' sparse-function
OPTS.use_mutils = 1;
    % 1 --> use MUTILS' "sparse_create"
OPTS.nthreads   = COMM.nthreads;
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


% =========================================================================
% PROBLEM SIZE, ALLOCATE SOLUTION VECTORS
% =========================================================================
nP         = MESH.nVnod;
nmg        = MESH.nmg;
nmesh      = length(MESH.EL2NOD);


% =========================================================================
% PARALLEL STUFF
% =========================================================================
nsd        = COMM.nsd;
if nsd>1
    NB        = COMM.NB;
    sumSDB_fh = COMM.sum_SDB;
    dot_fh    = COMM.dot;
    myPnodSDB = cell(1,COMM.nNB);
    for iNB=1:COMM.nNB
        mynod_SDB      = COMM.mynod_SDB{1}{iNB};
        myPnodSDB{iNB} = mynod_SDB( mynod_SDB<=nP );
    end
    unique_Pdofs = COMM.unique_nodes{1}( COMM.unique_nodes{1}<=nP );
    nP_D         = COMM.sum_all( length(unique_Pdofs) );
end
% compare_serial_parallel_3d(MESH,COMM,VAR.Dens,'Dens',0);
% compare_serial_parallel_3d(MESH,COMM,VAR.Visc,'Visc',1);

                                                                            tic
% =========================================================================
% ASSEMBLY OF GLOBAL MATRICES
% =========================================================================
switch SETTINGS.jacobian
    case 'standard'
        [K1,Fb,GG,MM_mu,Fpdil] = assembly_stokes_sph_3d...
            (VAR,MESH,PHYSICS,NUMSCALE,OPTS);                              Uinfo.Tasmble=toc;tic
        fprintf(fidl,' Assembly of %6i elements in %3i blocks   = %4.1f sec\n',...
            MESH.nel,ceil(MESH.nel/OPTS.nelblk),Uinfo.Tasmble);
    case 'double'
        [K1,Fb,GG,MM_mu,Fpdil] = assembly_stokes_sph_3d_double_jacobian...
            (VAR,MESH,PHYSICS,NUMSCALE,OPTS);                              Uinfo.Tasmble=toc;tic
        fprintf(fidl,' Assembly of %6i elements in %3i blocks   = %4.1f sec\n',...
            MESH.nel,ceil(MESH.nel/OPTS.nelblk),Uinfo.Tasmble);
end
% % [K1,Fb,GG,MM_mu,Fpdil] = assembly_stokes_sph_3d...
% %     (VAR,MESH,PHYSICS,NUMSCALE,OPTS);                                       Uinfo.Tasmble=toc;
% % fprintf(fidl,' Assembly of %6i elements in %3i blocks   = %5.1f sec\n',MESH.nel,ceil(MESH.nel/OPTS.nelblk),Uinfo.Tasmble);
% % 																			tic
% compare_serial_parallel_3d(MESH,COMM,diag(MM_mu),'diagMM',0);
% compare_serial_parallel_3d(MESH,COMM,Fb         ,'Fb',0);

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
    is_tril = full(all(K1(1,2:end)==0));
    K1      = K1 + tril(K1,-1)'; % Create upper half of symmetric K matrix
    K1      = RR*K1*RR';
    if is_tril
        K1 = tril(K1);
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
    is_tril = full(all(K1(1,2:end)==0));
    K1      = K1 + tril(K1,-1)'; % Create upper half of symmetric K matrix
    K1      = RR'*K1*RR;
    if is_tril
        K1 = tril(K1);
    end
    GG = RR' * GG;
    Fb = RR' * Fb;
end
                                                                            Uinfo.Trotate=toc;
fprintf(fidl,' Rotation of stiffness matrix                = %4.1f sec\n',Uinfo.Trotate);

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
rmPnullspace = 1; % IN A SPHERE THE DOMAIN IS CLOSED
% if nsd>1; rmPnullspace=COMM.minLabs(rmPnullspace); end
% if rmPnullspace
%     FigNo = 100 + FigNo;
    fprintf(fidl,' BCs DO NOT allow for unconstrained in- or outflow. Will remove nullspace.\n');
% else
%     fprintf(fidl,' BCs DO allow for unconstrained in- or outflow. Will NOT remove nullspace.\n');
% end

% =========================================================================
% BOUNDARY CONDITIONS
% =========================================================================
nU          = 3*MESH.nnod;
ifix        = UBC.ifix; % list of constrained (fixed) velocity degrees of freedom
Ufix        = UBC.vfix; % list of values at fixed velocity degrees of freedom
ifree       = 1:nU;
ifree(ifix) = [];       % list of free velocity degrees of freedom
U(ifix)     = Ufix;     % prescribed velocity values
U0          = U;

% Subtract boundary condition terms from r.h.s.:
% (if K1 would contain both upper and lower triangular parts it would 
% simply be: bc_rhs = K1(:,ifix)*Ufix(:);)
KD           = full(diag(K1));
bc_rhs       = K1(:,ifix)*Ufix(:);           % multiply lower part of K1
bc_rhs       = bc_rhs + K1(ifix,:)'*Ufix(:); % multiply "upper" part of K1
bc_rhs(ifix) = bc_rhs(ifix) - KD(ifix).*Ufix(:); % subtract diagonal (counted twice)

% KK1     = K1 + tril(K1,-1)';
% bc_rhs2 = KK1(:,ifix)*Ufix(:);

% =========================================================================
% CALCULATE R.H.S VECTOR
% =========================================================================
Rhs = Fb + GG*P - bc_rhs;
    % buoyancy forces + forces resulting from pressure gradients +
    % velocity boundary conditions
% compare_serial_parallel_3d(MESH,COMM,Rhs,'Rhs',0);
% compare_serial_parallel_3d(MESH,COMM,diag(K1),'diagK1',0);

Rhs = Rhs(ifree);
K1  = K1(ifree,ifree);
                                                                            tic
% =========================================================================
% PREPARATION OF STIFFNESS MATRICES ON ALL MULTIGRID LEVELS
% =========================================================================
% Restrict stiffness matrix to all multigrid levels
% Note: KK will be a cell array with "nmg" elements; each element
%       contains the stiffness matrix on the respective multigrid level.
[KK,Ic2f_U] = restrict_stiffmat(MESH,K1,ifree,fidl);                     Uinfo.Trestrict=toc;tic
fprintf(fidl,' Restriction of stiffness matrix             = %4.1f sec\n',Uinfo.Trestrict);

% If MUTILS is available: convert the sparse matrix to the more efficient
% format (see help sparse_convert)
if use_mutils && exist('sparse_convert','file')
    opts_spc.symmetric = 1;
    opts_spc.nthreads  = COMM.nthreads;
%     if nsd==1
%         opts_spc.nthreads = 8;
%         setenv('OMP_NUM_THREADS', num2str(opts_spc.nthreads));
%     end
    K1 = sparse_convert(K1,opts_spc);
    fprintf(fidl,' MUTILS sparse_convert on K                  = %4.1f sec\n',toc);tic
end

LcSD=[];permLcSD=[]; %#ok<NASGU>
if nsd>1
    % =====================================================================
    % Update communication data to account for the reduced size of system
    % of equations to be solved (ifree instead of nU). The system of 
    % equations is reduced on all multigrid levels.
    % =====================================================================
    COMM = include_BCs(COMM,MESH,ifix,ifree,nmg); % *SUBFUNCTION*
                                                                            Uinfo.Tbc2comm=toc;
    fprintf(fidl,' Modifying COMM structure for including BCs  = %4.1f sec\n',Uinfo.Tbc2comm);tic

    % =====================================================================
    % CHOLESKY FACTORIZATION OF STIFFNESS MATRIX OON BASE-LEVEL
    % =====================================================================
    % Merge the stiffness matrices of all subdomains on the BASE level and
    % perform a Cholesky factorization of stiffness matrix
%     [KbD,permKbD,Uinfo] = merge_KbD_v1(KK{nmesh},COMM,coarse_solver,fidl); % *SUBFUNCTION*
    [KbD,permKbD,Uinfo] = merge_KbD_v2(KK{nmesh},COMM,coarse_solver,fidl); % *SUBFUNCTION*
    
    % =====================================================================
    % This is required in parallel to correctly restrict the SD coarse mesh
    % residual to the BASE mesh (coarsest MG level)
    % =====================================================================
    if nmesh>nmg
        error(' This block needs to be verified');
        nUdof_c  = 3*MESH.nnod_mg(nmg);
        iUfree_c = ifree(ifree<=nUdof_c);
        i0       = setdiff(iUfree_c,COMM.unique_Udof{nmg});
        Ic2f_U{nmg}(i0,:) = 0;
    end
else
    switch coarse_solver
        case 'chol'
            % =====================================================================
            % CHOLESKY FACTORIZATION OF STIFFNESS MATRIX ON BASE-LEVEL
            % =====================================================================
            [KbD,permKbD] = cholesky_factorization(KK{nmesh});              Uinfo.t_chol=toc;
                                                                            tmp=whos('KbD');Uinfo.MB_KbD=tmp.bytes/1048576;
            fprintf(fidl,' Factorization of coarse stiffness matrix    = %4.1f sec (%1iMB)\n',...
                    Uinfo.t_chol,round(Uinfo.MB_KbD));
        case 'ichol'
            error(' to be coded');
        otherwise
            KbD     = KK{nmesh};
            permKbD = [];
    end
end


% =========================================================================
% PREPARATION OF PRESSURE PRECONDITIONER
% =========================================================================
if pc_Pat==1
    DMM = diag(MM_mu); % diagonal of MM_mu
    if nsd>1; DMM=sumSDB_fh(DMM,NB,myPnodSDB); end
    DiMM = 1./DMM;  % inverse diagonal of MM_mu
    clear MM_mu
elseif ismember(pc_Pat,[2 3])
    MM_mu = MM_mu + tril(MM_mu,-1)';
    LMM   = sum(MM_mu,2); % lumped MM_mu
    if nsd>1; LMM=sumSDB_fh(LMM,NB,myPnodSDB); end
    LiMM = 1./LMM;  % inverse lumped MM_mu
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

% compare_serial_parallel_3d(MESH,COMM,Fb,'Fb',0);
% compare_serial_parallel_3d(MESH,COMM,bc_rhs,'bc_rhs',0);
% compare_serial_parallel_3d(MESH,COMM,GG*P,'Fp',1);

if FigNo
    nsub=3;
    figure(FigNo);clf;
    set(gcf,'Renderer','zbuffer');
end

% =========================================================================
% CALCULATE INITIAL RESIDUAL VECTOR
% =========================================================================
OPTS_U = struct('useMG',useMG,...
                'itmax',itmax_U,...
                'tol'  ,rtol_U*rtol_restart,...
                'is_abs_tol',0,...
                'coarse_solver',coarse_solver);
% OPTS_U.useMG = 1;
% OPTS_U.itmax = 100;
% OPTS_U.FigNo = FigNo;
% save('Stokes_open_03','K1','KK','Rhs','U','ifree','Ic2f_U','nmg','COMM','OPTS_U');
[U(ifree),nit,time,rU_rms_vec] = mgpcg_solver_p...
    (K1,KK,KbD,permKbD,Rhs,U(ifree),Ic2f_U,nmg,COMM,OPTS_U);                nitU=nit;time_U=time;
OPTS_U.tol = rtol_restart*rtol_U;

% r = -rhs_P - S p + Fpdil
%   = -G'*Kinv*Fb - G'*Kinv*G*p + Fpdil , but Fp=G*p and U=Kinv*F (see above)
% 	= -G'*u_ib - G'*Kinv*Fp + Fpdil     , F=Fb+Fp
% 	= -G'*U + Fpdil
r     = -GG'*U + Fpdil; % initial pressure residual vector
% compare_serial_parallel_3d(MESH,COMM,U,'U',0);
% compare_serial_parallel_3d(MESH,COMM,r,'r',1);

if nsd>1; r = sumSDB_fh(r,NB,myPnodSDB); end
if rmPnullspace; r = rm0space_p(r); end
norm_rP = normdf_p(r); % norm of residual vector (shown in convergence plot)
tol_Pat = rtol_Pat*norm_rP;

% Solve again but this time with absolute tolerance depending on tol_Pat
OPTS_U.tol        = rtol_U*rtol_restart*norm_rP;
OPTS_U.is_abs_tol = 1;
[U(ifree),nit,time,rms_vec] = mgpcg_solver_p...
    (K1,KK,KbD,permKbD,Rhs,U(ifree),Ic2f_U,nmg,COMM,OPTS_U);                nitU=nitU+nit;time_U=time_U+time;rU_rms_vec(end+1:end+nit+1)=rms_vec;
OPTS_U.tol        = rtol_U*tol_Pat;

% Plot convergence
% ================
if FigNo
    figure(FigNo);
    ah1=subplot(1,nsub,1); plot(0,log10(norm_rP),'ro');
    if ~isempty(YLimit); set(gca,'YLim',YLimit); end
    hold all
    set(gca,'FontSize',fs); grid on
    xlabel('P iteration'); ylabel('Norm of pressure residual vector');
    line([0 itmax_Pat],log10([tol_Pat tol_Pat]),...
          'Linestyle','-','Color','r','Linewidth',2);
    set(gca,'XLim',[0 itmax_Pat]); hold all

    ah2=subplot(1,nsub,2); plot(log10(rU_rms_vec),'r.-');
    if ~isempty(YLimit); set(gca,'YLim',YLimit); end
    hold on
    set(gca,'FontSize',fs); grid on
    xlabel('Z iteration'); ylabel('Norm of velocity residual vector');
%     line([0 itmax_U],log10([tol_U tol_U]),...
%           'Linestyle','--','Color','k','Linewidth',2);
    set(gca,'XLim',[0 itmax_U]); hold all
    
    ah3=subplot(1,nsub,3); plot(0,0,'w.');
    if ~isempty(YLimit); set(gca,'YLim',YLimit); end
    hold on
    set(gca,'FontSize',fs); grid on
    xlabel('U iteration'); ylabel('Norm of velocity residual vector');
%     line([0 itmax_U],log10([tol_U tol_U]),...
%           'Linestyle','--','Color','k','Linewidth',2);
    set(gca,'XLim',[0 itmax_U]); hold all
    
%     saveas(gcf,[SETTINGS.outdir '/' COMM.prefix '_StokesFlowSolverConvergence'],'png');
end

% % Show profiling data in terminal
% % ===============================
% fprintf(' Absolute tolerances to calculate initial rP = %0.2e\n',tol_U0);
% fprintf(' Relative/absolute tolerances for pressure   = %0.2e / %0.2e\n',rtol_Pat,tol_Pat);
% fprintf(' Relative/absolute tolerances for velocity   = %0.2e / %0.2e\n',rtol_U,tol_U);

dU = max(abs(U0-U)); clear U0
if nsd>1; dU = COMM.maxLabs(dU); end

display_profiling_data(0);

nitZ       = 0;
time_Z     = 0;
restart_CG = 1;
convergence_problem = 0;
                                                                            time_Pat=0;
% =========================================================================
%  BEGIN OF PRESSURE ITERATIONS
% =========================================================================
for itPat = 1:itmax_Pat
                                                                            tic
    % Pressure preconditioning (approximate error of current P)
    % =========================================================
    d = precondition_P; % *NESTED FUNCTION*

    % get new search direction q
    % ==========================
    if restart_CG
        % Initialize Conjugate Gradients
        % ==============================
        q          = d; % define FIRST search direction q
        restart_CG = 0;
        norm_rP0   = norm_rP;
        if FigNo
            plot([0 itmax_Pat],log10([rtol_restart*norm_rP0 rtol_restart*norm_rP0]),...
                'b--','Linewidth',2,'Parent',ah1);
        end
        OPTS_Z = struct('useMG',useMG,...
                        'itmax',itmax_U,...
                        'tol'  ,rtol_U*max(rtol_restart*norm_rP0,tol_Pat),...
                        'is_abs_tol',1,...
                        'coarse_solver',coarse_solver);
    else
        rdlast = rd; % save old rd (for calculating beta)
        % Make new search direction q S-orthogonal to all previous q's
        beta   = dot_p(r-rlast,d)/rdlast;  % Polak-Ribiere version 1
        q      = d + beta*q; % calculate NEW search direction
    end
    
    if rmPnullspace; q = rm0space_p(q); end % *SUBFUNCTION*
    rd = dot_p(r,d); % numerator for calculating alpha

    % Perform the S times q multiplication
    % ====================================
    % S cannot be calculated explicitly, since Kinv cannot be formed
    % Sq = S * q = (G' * Kinv * G) * q
    % Hence, the muliplicatipon is done in 3 steps:
    % (1) y   = G*q
    % (2) K z = y (Multigrid-preconditioned Conjugate Gradient solver)
    % (3) Sq  = G'*z
    y  = GG*q;                                                              time_Pat=time_Pat+toc;
    z  = zeros(nU,1);
    [z(ifree),nit,time,rZ_rms_vec] = mgpcg_solver_p...
        (K1,KK,KbD,permKbD,y(ifree),z(ifree),Ic2f_U,nmg,COMM,OPTS_Z);       time_Z=time_Z+time;nitZ=nitZ+nit;
    if FigNo
        plot([0 itmax_U],log10([OPTS_Z.tol OPTS_Z.tol]),...
            'b--','Linewidth',2,'Parent',ah2);
    end
                                                                            tic
    Sq = GG'*z;
    if nsd>1; Sq = sumSDB_fh(Sq,NB,myPnodSDB); end
    
    qSq   = dot_p(q,Sq); % denominator in calculating alpha
    rlast = r;      % needed for Polak-Ribiere version 1
    alpha = rd/qSq; % step size in direction q
    
    % Update solution and residual
    % ============================
    P = P + alpha*q;  % update pressure solution
    U = U + alpha*z;  % update velocity by accumulating alpha*z
    r = r - alpha*Sq; % update residual vector
    if rmPnullspace; r = rm0space_p(r); end % remove nullspace from both q AND Sq
    norm_rP_last = norm_rP;
    norm_rP      = normdf_p(r);  % norm of residual vector (shown in convergence plot)
    norm_rP_min  = min(norm_rP_last,norm_rP);
                                                                            time_Pat=time_Pat+toc;
    % Show profiling data, update convergence plot
    % ============================================
    if FigNo
        plot(itPat,log10(norm_rP),'k.','Parent',ah1);
        plot(log10(rZ_rms_vec),'k.-','Parent',ah2); drawnow;
%         saveas(gcf,[SETTINGS.outdir '/' COMM.prefix '_StokesFlowSolverConvergence'],'png');
    end
    
    dU = max(abs(alpha*z));
    if nsd>1; dU = COMM.maxLabs(dU); end
    display_profiling_data(itPat);
    
    % Check if solution converged
    % ===========================
    if norm_rP<tol_Pat && itPat>=itmin_Pat
        break
    end
    
    % Check if restart is required
    % ============================
    if norm_rP < rtol_restart*norm_rP0 || norm_rP>10*norm_rP0
        % Calculate exact velocity for current pressure using accumulated
        %  "U" as initial guess; calculate true divergence
        % ===================================================================
        Rhs = Fb + GG*P - bc_rhs;
        % compare_serial_parallel_3d(MESH,COMM,Rhs,'Rhs',0);
        [U(ifree),nit,time,rU_rms_vec] = mgpcg_solver_p...
            (K1,KK,KbD,permKbD,Rhs(ifree),U(ifree),Ic2f_U,nmg,COMM,OPTS_U); time_U=time_U+time;nitU=nitU+nit;tic
        r  = -GG'*U + Fpdil;
        if nsd>1; r = sumSDB_fh(r,NB,myPnodSDB); end
        % compare_serial_parallel_3d(MESH,COMM,r,'r',0);
        
        norm_rP_true = normdf_p(r);  % norm of residual vector (shown in convergence plot)
        if norm_rP_true>1.1*norm_rP
            rtol_U = rtol_U * 0.75*(norm_rP/norm_rP_true);
        end
                                                                            time_Pat=time_Pat+toc;
        if norm_rP>10*norm_rP0 || norm_rP>1.2*norm_rP_min
            if convergence_problem
                break
            end
            rtol_U = 0.1*rtol_U; % we had a convergence problem; increase accuracy of velocity solver
            convergence_problem = 1;
        end
        
        if FigNo
            plot(itPat,log10(norm_rP_true),'bo','Parent',ah1);
            plot(log10(rU_rms_vec),'b.-','Parent',ah3); drawnow;
%             saveas(gcf,[SETTINGS.outdir '/' COMM.prefix '_StokesFlowSolverConvergence'],'png');
        end
        restart_CG = 1;
        % compare_serial_parallel_3d(MESH,COMM,U,'U',1);
    end

end % END OF PRESSURE ITERATIONS

% % % Store pressure residual
% % VAR.rP = r;

% compare_serial_parallel_3d(MESH,COMM,r,'r',0);
% compare_serial_parallel_3d(MESH,COMM,P,'P',0);
% compare_serial_parallel_3d(MESH,COMM,U,'U0',0);

% Solve for U one last time
Rhs = Fb + GG*P - bc_rhs;
[U(ifree),nit,time,rU_rms_vec] = mgpcg_solver_p...
    (K1,KK,KbD,permKbD,Rhs(ifree),U(ifree),Ic2f_U,nmg,COMM,OPTS_U);         time_U=time_U+time;nitU=nitU+nit;tic
if FigNo
    plot(log10(rU_rms_vec),'r.-','Parent',ah3); drawnow;
    saveas(gcf,[SETTINGS.outdir '/' COMM.prefix '_StokesFlowSolverConvergence'],'png');
end

% % % Store true divergence of last flow field
% % VAR.divU = -GG'*U + Fpdil;

% if norm_rP>tol_Pat
%     fprintf('\n WARNING: Pressure iterations did not converge within %3i iterations: tol=%0.2e norm_rP=%0.2e\n',...
%             itPat,tol_Pat,norm_rP);
% end

fprintf(fidl,'\n Time spent on pressure problem               = %6.1f sec\n',time_Pat);
fprintf(fidl,' Time spent on velocity problem               = %6.1f sec\n',time_U);
fprintf(fidl,' Total number of Patera iterations            = %5i\n',itPat);
fprintf(fidl,' Total number of Z iterations                 = %5i\n',nitZ);
fprintf(fidl,' Total number of velocity iterations          = %5i\n',nitU);
fprintf(fidl,' Total time for solving flow problem          = %6.1f sec\n',etime(clock,clock0));
fprintf(fidl,'\n');

if FigNo
    figure(FigNo);
    saveas(gcf,[SETTINGS.outdir '/' COMM.prefix '_StokesFlowSolverConvergence'],'png');
end

% compare_serial_parallel_3d(MESH,COMM,U,'Rhs',0);
% compare_serial_parallel_3d(MESH,COMM,U,'U',1);

% Check solution
if any(isnan(U)) || any(isnan(P))
	error('Solution containes NaNs !!!!');
end

if strcmp(SETTINGS.plate_model,'WJM')
    VAR.Ur = U(1:3:end); % Radial component
    VAR.Ue = U(2:3:end); % East component
    VAR.Un = U(3:3:end); % North component
    U_xyz  = RR'*U;      % rotate velocity solution back to xyz-coordinates
elseif any(strcmp(SETTINGS.plate_model,{'GPlates','Debug'}))
    VAR.Uth = U(1:3:end); % theta component
    VAR.Uph = U(2:3:end); % phi component
    VAR.Ur  = U(3:3:end); % radial component
    U_xyz   = RR*U;       % rotate velocity solution back to xyz-coordinates
end
VAR.Ux = U_xyz(1:3:end); % x-velocity solution
VAR.Uy = U_xyz(2:3:end); % y-velocity solution
VAR.Uz = U_xyz(3:3:end); % z-velocity solution
VAR.P  = P;          % pressure solution

% compare_serial_parallel_3d(MESH,COMM,U,'U',0);
% compare_serial_parallel_3d(MESH,COMM,P,'P',1);

% Save profiling data
Uinfo.tP     = time_Pat;
Uinfo.tU     = time_U;
byte2MB      = 1/1048576;
KK1          = KK{1}; %#ok<NASGU>
tmp          = whos('KK1'); clear KK1
Uinfo.memKK  = tmp.bytes * byte2MB;
tmp          = whos('KbD');
Uinfo.memLL  = tmp.bytes * byte2MB;
tmp          = whos('GG');
Uinfo.memGG  = tmp.bytes * byte2MB;
Uinfo.nel    = MESH.nel;
Uinfo.nP     = MESH.nVnod;
Uinfo.nU     = nU;
Uinfo.nfree  = length(ifree);
Uinfo.nitP   = itPat;
Uinfo.norm_rP = norm_rP;
Uinfo.nitU   = nitU;

% =========================================================================
%                            NESTED FUNCTIONS
% =========================================================================

function d = precondition_P
    % Get approximation d to the unknown error of the current
    % pressure solution P (i.e. estimate approx error d from residual r)
    switch pc_Pat
        case 0 % No preconditioning (i.e. CG algorithm)
            d = r;
        case 1 % Scaling by diagonal of (1/viscosity)-scaled mass matrix
            d = r.*DiMM;
        case 2 % Scaling by lumped (1/viscosity)-scaled mass matrix
            d = r.*LiMM;
        case 3 % Jacobi iterations on (1/viscosity)-scaled mass matrix
            d  = r .* LiMM;
            for is=1:length(wght_MM)
                Md = MM_mu*d;
                if nsd>1; Md = sumSDB_fh(Md,NB,myPnodSDB); end
                d  = d + wght_MM(is) * LiMM .* (r - Md);
            end
    end
end % END OF NESTED FUNCTION precondition_P

% =========================================================================

function display_profiling_data(flag)
if flag==0
    fprintf(fidl,'               Patera algorithm              \n');
    fprintf(fidl,'  itP |  time   |    rP    |    rU    |    dU    | nitU \n');
    fprintf(fidl,' %4i | %7.1f | %0.2e | %0.2e | %0.2e | %4i \n',...
            0,etime(clock,clock0),norm_rP,rU_rms_vec(end),dU,nit);
else
    fprintf(fidl,' %4i | %7.1f | %0.2e | %0.2e | %0.2e | %4i \n',...
            itPat,etime(clock,clock0),norm_rP,rU_rms_vec(end),dU,nit);
end
end % END OF NESTED FUNCTION display_profiling_data

% =========================================================================

function rms = normdf_p(a)
    if nsd>1
        rms = COMM.normdf(a,unique_Pdofs);
%         rms = dot_fh(a,a,unique_Pdofs);
%         rms = sqrt(rms/nP_D);
    else
        rms = normdf(a);
    end
end % END OF NESTED FUNCTION normdf_p

% =========================================================================

function ab = dot_p(a,b)
    if nsd>1
        ab = dot_fh(a,b,unique_Pdofs);
    else
        ab = dot(a,b);
    end
end % END OF NESTED FUNCTION dot_p

% =========================================================================

function a = rm0space_p(a)
    % Function takes out the mean of vector a
    % (thereby removes the constant pressure mode, i.e. the nullspace)
    if nsd>1
        sum_a  = sum( a(unique_Pdofs) );
        mean_a = COMM.sum_all( sum_a ) / nP_D;
    else
        % nDpbc     = length(a2); use to check this variable in parallel
        mean_a = mean(a);
    end
    a = a - mean_a;
end % END OF NESTED FUNCTION rm0space_p

end % END OF FUNCTION flow3d_patera_mg_sph_p

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function [KbD,permKbD,Uinfo] = merge_KbD_v2(Kb,COMM,coarse_solver,fidl)

myid = COMM.myid;
nsd  = COMM.nsd;

% Get sparse pattern of BASE stiffness matrix Kb in this subdomain
% "Kbi" and "Kbj" are the locations of the nonzero values "Kbv"
% within the sparse matrix Kb of coarsest mesh (BASE)
[Kbi,Kbj,Kbv] = find(tril(Kb)); % only lower half of Kc is needed
Kbi           = int32(Kbi);
Kbj           = int32(Kbj);
nSDdofs       = length(Kbi);    % number of non-zero entries
myUdof_SD2D   = COMM.Udof_SD2D_c{myid};

% Sum up the number of non-zeros of all subdomains to allocate memory
nDdofs = COMM.sum_all(nSDdofs);
if myid==1
    % Allocate memory for triple format indices of whole-domain coarse KK
    KbiD_all = zeros(nDdofs,1,'int32');
    KbjD_all = zeros(nDdofs,1,'int32');
    KbvD_all = zeros(nDdofs,1);

    KbiD     = myUdof_SD2D(Kbi); clear Kbi
    KbjD     = myUdof_SD2D(Kbj); clear Kbj
    ii       = 1;
    n        = length(KbiD); % same as length(jDdofs)
    KbiD_all(ii:ii+n-1) = KbiD(:); clear KbiD
    KbjD_all(ii:ii+n-1) = KbjD(:); clear KbjD
    KbvD_all(ii:ii+n-1) = Kbv(:);  clear Kbv
    ii                  = ii + n;
                                                            MB_Kb=0;tic
end

% Now all labs send the triplet of their subdomain Kb to SD 1. The pointer
% myUdof_SD2D is used to translate the subdomain dofs into domain dofs.
for isd=2:nsd
    labBarrier
    tag = 10000*isd+1;
    if myid==1
        % I receive domain dofs (rows) from SD "isd"
        KbiD = COMM.recv_1NB(isd,tag);
    elseif myid==isd
        % I send my domain dofs (rows) to SD "1"
        COMM.send_1NB(1,myUdof_SD2D(Kbi),tag); clear Kbi
    end
    labBarrier
    tag = 10000*isd+2;
    if myid==1
        % I receive domain dofs (columns) from SD "isd"
        KbjD = COMM.recv_1NB(isd,tag);
    elseif myid==isd
        % I send my domain dofs (columns) to SD "1"
        COMM.send_1NB(1,myUdof_SD2D(Kbj),tag); clear Kbj
    end
    labBarrier
    tag = 10000*isd+3;
    if myid==1
        % I receive domain dofs (columns) from SD "isd"
        KbvD = COMM.recv_1NB(isd,tag);
    elseif myid==isd
        % I send my domain dofs (columns) to SD "1"
        COMM.send_1NB(1,Kbv,tag); clear Kbv
    end

    if myid==1
                                                                            tmp=whos('KbvD');MB_Kb=MB_Kb+tmp.bytes/1048576;
                                                                            tmp=whos('KbiD');MB_Kb=MB_Kb+2*tmp.bytes/1048576;
        % Construct the triplet of the whole domain coarse stiffness matrix
        if all(KbiD==KbjD)
            error('KbiD and KbjD are identical! Communication screwed up?')
        end
        n                   = length(KbiD); % same as length(jDdofs)
        KbiD_all(ii:ii+n-1) = KbiD(:); clear KbiD
        KbjD_all(ii:ii+n-1) = KbjD(:); clear KbjD
        KbvD_all(ii:ii+n-1) = KbvD(:); clear KbvD
        ii                  = ii + n;
    end
end
                                                                            Uinfo.t_broadcast_Kb=toc;
if myid==1
                                                                            Uinfo.MB_Kb=MB_Kb;
    fprintf(fidl,' Received coarse subdomain stiffness matrices  = %4.1f sec (%1iMB)\n',...
                    Uinfo.t_broadcast_Kb,round(MB_Kb));

    % Create the DOMAIN BASE stiffness matrix KbD using the above triplet
    KbD = sparse3(KbiD_all,KbjD_all,KbvD_all);                              tmp=whos('KbD');Uinfo.MB_KbD=tmp.bytes/1048576;Uinfo.t_KbD=toc;
    % Get rid of arrays that are no longer required
    clear KbiD_all KbjD_all KbvD_all
    fprintf(fidl,' Assembly of global coarse stiffness matrix    = %4.1f sec (%1iMB)\n',...
                     Uinfo.t_KbD,round(Uinfo.MB_KbD)); 

    if strcmp(coarse_solver,'chol')
    % Perform Cholesky factorization of DOMAIN stiffness matrix on coarses mesh
    % =========================================================================
    [KbD,permKbD] = cholesky_factorization(KbD);                            Uinfo.t_chol=toc;tmp=whos('KbD');Uinfo.MB_KbD=tmp.bytes/1048576;
    fprintf(fidl,' Factorization of coarse stiffness matrix      = %4.1f sec (%1iMB)\n',...
                     Uinfo.t_chol,round(Uinfo.MB_KbD));
    else
    permKbD = [];
    end
else
    KbD     = [];
    permKbD = [];
end

end % END OF SUBFUNCTION merge_KbD_v2

% ###################################################################

function COMM = include_BCs(COMM,MESH,ifix,ifree,nmg)

myid  = COMM.myid;
nsd   = COMM.nsd;
nNB   = length(COMM.NB);
nmesh = length(MESH.EL2NOD);
COMM.myUdofSDB   = cell(1,nmg);
COMM.unique_Udof = cell(1,nmg);
for img=1:nmg
    nU     = 3*MESH.nnod_mg(img);
    ifree = ifree(ifree<=nU);
    nfree = length(ifree);
    ifix  = ifix(ifix<=nU);
    if nU-(nfree+length(ifix))~=0
        error(' ifix or ifree are not correct');
    end

    % construct pointer from "all" (free+fixed) dofs to "free" dofs
    all2free_Udofs         = zeros(1,nU);
    all2free_Udofs(ifree) = 1:nfree;

    % use the pointer to extract and re-number the free dofs
    mynod_SDB_MGi = COMM.mynod_SDB{img};
    myUdofSDB     = cell(1,nNB);
    for iNB=1:nNB
        if COMM.NB(iNB)
            % Extract the SD boundary dofs shared with SD "NB(iNB)" but without
            % changing their order!
            myUdofSDB_MGi = nod2dof(mynod_SDB_MGi{iNB},3);
            tmp           = myUdofSDB_MGi(~ismember(myUdofSDB_MGi,ifix));
            % Use the above pointer to re-number the SD boundary dofs
            % The dofs now correspond to the reduced system of equations:
            % KK(free,free) * u(free) = rhs(free)
            myUdofSDB{iNB} = all2free_Udofs( tmp );
        end
    end
    COMM.myUdofSDB{img} = myUdofSDB;

    % Same re-numbering for the list of uniq dofs on finest MG level
    uniqUdofs             = nod2dof(COMM.unique_nodes{img},3);
    uniqUdofs_free        = uniqUdofs(~ismember(uniqUdofs,ifix));
    uniqUdofs_free        = all2free_Udofs( uniqUdofs_free );
    COMM.unique_Udof{img} = uniqUdofs_free;
end

% Calculate the number of free velocity dofs in whole domain
% (needed for norm calculation below)
COMM.nUD = zeros(1,nmesh);
for img=1:nmg
    COMM.nUD(img) = COMM.sum_all( length(COMM.unique_Udof{img}) );
end
if nmesh>nmg
    COMM.nUD(nmesh) = 3*MESH.nnod(nmesh);
end

% Have to also change the pointer SD --> domain for the velocity dofs
% on the coarsest mesh. This is only required, if the coaresest mesh is
% the a geometic multigrid mesh (i.e. parent to the next finer mesh).
% If, however, the coarest mesh is an extra-coarse BASE mesh (i.e.
% nmg<nmesh) this is not necessary.
if nmesh==nmg
    % Pointer SD dofs --> domain dofs
    myUdof_SD2D_c = nod2dof(COMM.nod_SD2Dc{myid},3);

    % get a unique list of the free domain dofs on the coarse mesh
    ifree_D = COMM.unique(myUdof_SD2D_c(ifree))';

    % construct pointer all-to-free domain dofs on the coarsest MG level
    all2free_Udofs_D           = zeros(1,3*MESH.nnod_D);
    all2free_Udofs_D(ifree_D) = 1:length(ifree_D);

    myUdof_SD2D_c = all2free_Udofs_D( myUdof_SD2D_c(ifree) );
        
% if numlabs>1
%     filename = ['Workspace2_Lab' num2str(labindex)];
%     save(filename);
%     stop
% else
%     figure(40);clf;
%     for isd=1:nsd
%         filename = ['Workspace2_Lab' num2str(isd)];
%         load(filename);
%         showFreeDofs_SD2D(40,MESH{nmg}.gX,MESH{nmg}.EL2NOD,myUdof_SD2D_c);
%     end
% end
        
% if numlabs==1; return; end

    Udof_SD2D_c = cell(1,nsd);
    for isd=1:nsd
        if myid==isd
            Udof_SD2D_c{isd} = myUdof_SD2D_c;
            labBroadcast( myid , myUdof_SD2D_c );
        else
            Udof_SD2D_c{isd} = labBroadcast( isd );
        end
    end
    COMM.Udof_SD2D_c = Udof_SD2D_c;
end

end % END OF SUBFUNCTION include_BCs