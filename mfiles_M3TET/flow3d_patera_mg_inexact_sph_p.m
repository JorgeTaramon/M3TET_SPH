function [VAR,MESH,Uinfo] = flow3d_patera_mg_inexact_sph_p(VAR,MESH,COMM,UBC,SETTINGS,PHYSICS,NUMSCALE,dt)
% Usage: [VAR,MESH,Uinfo] = flow3d_patera_mg_inexact_sph_p(VAR,MESH,COMM,UBC,SETTINGS,PHYSICS,NUMSCALE,dt)
%
% Purpose: Stokes flow solver (calculate velocity and pressure solution)
% Input:
%   VAR      : [structure] : major variable fields, each is a vector
%   MESH     : [structure] : FE mesh parameters
%   COMM     : [structure] : inter-subdomain communication data
%   UBC      : [structure] : velocity boundary conditions
%   SETTINGS : [structure] : model parameters
%   PHYSICS  : [structure] : physical properties
%   NUMSCALE : [structure] : numerical scaling parameters
%   dt       : [double]    : time step
%
% Output:
%   VAR      : [structure] : structure containing all major variables
%   MESH     : [structure] : FE mesh parameters
%   Uinfo    : [structure] : profiling data
%
% Part of 3D convection code M3TET_MG, which is a simplified version of the
% parallel 3D mantle convection code M3TET (developed by
% J.Hasenclever & J.Phipps Morgan, 2007-2010)
% Email contact: joerg.hasenclever@zmaw.de
% For numerical methods see online Ph.D. thesis
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% Details on the Stokes flow solver:
%   Patera pressure iterations, serial version
%   Elements: 3D Taylor-Hood triangles (10 velocity nodes, 4 pressure nodes)
%   Viscous flow solution: Conjugate Gradient solver, preconditioned by a
%                          single V-cylce of geometric multigrid
%   Pressure solution: Preconditioned Conjugate Gradient (PCG) solver
%   Preconditioners: 0) none (algorithm reduces to CG)
%                    1) diagonal of (1/viscosity)-scaled mass matrix (DMM)
%                    2) d = inv(MM) r (testing only !!)
%                    3) Jacobi iterations on (1/viscosity)-scaled mass matrix
%   Null-space removed by taking out the mean of the pressure solution (if
%     domain is closed, i.e. no unconstrained in- or outflow)
%   Flexible CG used to orthogonalize q because of inexact preconditioner
%   Element assembly using standard method OR block-assembly (MILAMIN style) 
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
%
                                                                            clock0=clock;tic
% =========================================================================
% SETUP FOR STOKES FLOW SOLVER
% =========================================================================
FigNo    = 31;
    % number of figure showing convergence; ==0 --> no plot
fs       = 10;
    % font size in convergence plot
YLimit   = [-10 14];
    % y-axis limits in convergence plot
K_store = 'tril';
    % 'tril' --> store only lower triangular part of symmetric stiffness
    %            matrix (less memory usage but slower multiplications)
    % 'full'  --> store both halfes of stiffness matrix
% PART 1 - PRESSURE ITERATIONS
rtol_Pat   = 1e-5;
    % RELATIVE tolerance for pressure solution
    % (relative to norm of initial pressure residual)
PC_Pat     = 4;
    % 0 ==> no preconditioner (d=r)
    % 1 ==> d = r./diag(MM);  scaling by diagonal of mass matrix
    % 2 ==> d = r./LMM; scaling by lumped mass matrix; LMM=sum(MM,2)
    % 3 ==> d = d + (r - MM*d)./LMM; Jacobi iterations on MM
    % 4 ==> inexact Patera algorithm solving S d = r
nit_MM     = 12;
    % number of mass matrix iterations when using PC_Pat==3 or PC_Pat==4
if PC_Pat==4
    itmin_Pat = 1;   % minimum number of pressure iterations
    itmax_Pat = 15;  % maximum number of pressure iterations
else
    itmin_Pat = 3;   % minimum number of pressure iterations
    itmax_Pat = 120; % maximum number of pressure iterations
end

% PART 2 - PRECONDITIONING OF INEXACT PATERA ALGORITHM (PC_Pat==4)
OPTS_PC.rtol   = 2e-2;
    % tolerance of preconditioning (inexact) Patera relative to current 
    % norm of pressure residual
OPTS_PC.itmin  = 5;
    % minimum number of inexact Patera iterations
OPTS_PC.itmax  = 150;
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
dU_max     = 0.01;
    % Exit iterations if all nodal velocity components change by less than
    % this fraction
coarse_solver = 'chol'; % use "chol" !!!!!!!!!
    % Choose solver on coarsest multigrid level
    % 'chol'  --> Cholesky direct solver
    % 'ichol' --> incomplete Cholesky factorization
    % 'agmg'  --> Yvan Notay's algebraic multigrid solver AGMG
use_mutils = 1;
    % 1 --> use MUTILS tools for better preformance
    % 0 --> use Matlab's functions (slower)

% PART 4 - ELEMENT ASSEMBLY
OPTS.method    = 'opt';
    % std = standard element assembly (slow)
    % opt = optimized assembly (blocks of elements; MILAMIN method)
    %       see Dabrowski et al. 2008 (doi:10.1029/2007GC001719)
% when using MILAMIN method:
OPTS.nelblk    = 1000;
    % when using 'opt': number of elements assembled at once (block size)
    % Best speed depends on computer architecture (CPU cache)
    % try numbers between 200 and 1000+
OPTS.nel_asmbl = 1e8;
    % max number of elements for sparse command
OPTS.use_mutils = 1;
    % 1 --> use MUTILS' "sparse_create"
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
assembly_order = 1;
if ~isempty(strfind(get_hostname,'node')) && ...
   COMM.nsd>1 && COMM.maxLabs(MESH.nnod)>150000
    % Abakus has insufficient memory per CPU, which leads to swaping during
    % the element assembly. One way is to 
    assembly_order = 2-mod(COMM.myid,2);
end
for i=1:2
    if i==assembly_order
        [KK,Fb,GG,MM_mu,Fpdil,DensEl,ViscIP] = assembly_stokes_sph_3d...
            (VAR,MESH,PHYSICS,NUMSCALE,OPTS);                               Uinfo.Tasmble=toc;
    end
    labBarrier
end
fprintf(' Assembly of %6i elements in %3i blocks   = %5.1f sec\n',MESH.nel,ceil(MESH.nel/OPTS.nelblk),Uinfo.Tasmble);

KK = KK + tril(KK,-1)'; % Create upper half of symmetric K matrix
RR = make_rot_matrix(MESH,UBC.PlateInfo,SETTINGS.useNataf);
KK = RR*KK*RR'; %***WATCH OUT!!!! maybe RR should be RR' and vice versa
if strcmp(K_store,'tril')
    KK = tril(KK);
end
GG = RR*GG;
Fb = RR*Fb;
                                                                            tic
% compare_serial_parallel_3d(MESH,COMM,diag(MM_mu),'diagMM',0);
% compare_serial_parallel_3d(MESH,COMM,Fb         ,'Fb',0);


% =========================================================================
% FREE SURFACE TERMS
% =========================================================================
switch SETTINGS.top_surface
    case 'free'
        if dt>0
            % Free surface
            alpha    = SETTINGS.fs_alpha; % weighting begin/end of time step
            KKfs_zz  = free_surface_terms(MESH,PHYSICS,NUMSCALE,DensEl,alpha,dt);
            switch K_store
                case 'full'
                    KK = KK - KKfs_zz; 
                case 'tril'
                    KK = KK - tril(KKfs_zz);
                otherwise
                    error(' "K_store" must be either "tril or "full".');
            end
            clear KKfs_zz
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
nP          = length(P);
U0          = U;

switch K_store
    case 'full' % KK contains both upper and lower triangular part
        bc_rhs       = KK(:,ifix)*Ufix(:);

    case 'tril' % KK is only lower half
        bc_rhs       = KK(:,ifix)*Ufix(:);
            % multiply lower part of KK
        bc_rhs       = bc_rhs + KK(ifix,:)'*Ufix(:);
            % multiply "upper" part of KK
        KD           = full(diag(KK));
        bc_rhs(ifix) = bc_rhs(ifix) - KD(ifix).*Ufix(:);
            % subtract diagonal because it was counted twice
        
    otherwise
         error(' "K_store" must be either "tril or "full".');
end


% =========================================================================
% PARALLEL STUFF
% =========================================================================
nmg        = MESH.nmg;
nmesh      = length(MESH.EL2NOD);
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
rmPnullspace = 1; %is_domain_closed(MESH.PointID,UBC.ifix); %,MESH.GCOORD,MESH.EL2NOD{1});
if nsd>1; rmPnullspace=COMM.minLabs(rmPnullspace); end
if rmPnullspace
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
fprintf(' Restriction of stiffness matrix             = %4.1f sec\n',Uinfo.Trestrict);

% If MUTILS are available: convert the sparse matrix to more efficient
% format (see help sparse_convert)
if use_mutils
    addpaths_mutils();
    % if exist('sparse_convert','file')
    opts_spc.symmetric = 1;
%     if nsd==1
%         opts_spc.nthreads = 8;
%         setenv('OMP_NUM_THREADS', num2str(opts.nthreads));
%     end
    opts_spc.nthreads  = 1;
    K1 = KK{1};
    K1 = sparse_convert(K1,opts_spc);
else
    K1 = KK{1};
end
KK{1} = full(diag(KK{1}));

% if numlabs>1
%     save([COMM.prefix '_Workspace1']);
% else
%     load('Lab02x01_Workspace1');
% end

LcSD=[];permLcSD=[]; %#ok<NASGU>
if nsd>1
    
% if numlabs>1
    % =====================================================================
    % Update communication data to account for the reduced size of system
    % of equations to be solved (ifree instead of nU). The system of 
    % equations is reduced on all multigrid levels.
    % =====================================================================
    COMM = include_BCs(COMM,MESH,ifix,ifree,nmg); % *SUBFUNCTION*
                                                                            Uinfo.Tbc2comm=toc;
    fprintf(' Modifying COMM structure for including BCs  = %4.1f sec\n',Uinfo.Tbc2comm);tic

%     save([COMM.prefix '_Workspace2']);
%     stop
% else
%     load('Lab02x01_Workspace2');
% end

    % =====================================================================
    % CHOLESKY FACTORIZATION OF STIFFNESS MATRIX OON BASE-LEVEL
    % =====================================================================
    % Merge the stiffness matrices of all subdomains on the BASE level and
    % perform a Cholesky factorization of stiffness matrix
%     [KbD,permKbD,Uinfo] = merge_KbD_v1(KK{nmesh},COMM,coarse_solver); % *SUBFUNCTION*
    [KbD,permKbD,Uinfo] = merge_KbD_v2(KK{nmesh},COMM,coarse_solver); % *SUBFUNCTION*
    
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
    case 'ichol'
        % =====================================================================
        % INCOMPLETE CHOLESKY FACTORIZATION OF STIFFNESS MATRIX ON BASE-LEVEL
        % =====================================================================
        clear opts_ichol
%         opts_ichol.type     = 'nofill';
        opts_ichol.type     = 'ict';
        opts_ichol.droptol  = 1e-14 * max(max(abs(KK{nmesh})));
        opts_ichol.michol   = 'on';
        opts_ichol.diagcomp = max(sum(abs(KK{nmesh}),2)./diag(KK{nmesh}))-2;
        permKbD             = amd(KK{nmesh});
        KbD                 = cs_transpose(KK{nmesh});
        KbD                 = cs_symperm(KbD,permKbD);
        KbD                 = cs_transpose(KbD);
        KbD                 = ichol(KbD,opts_ichol);
%         tmp=whos('KbD');tmp.bytes/1024^2
    case 'chol'
        % =====================================================================
        % CHOLESKY FACTORIZATION OF STIFFNESS MATRIX ON BASE-LEVEL
        % =====================================================================
        [KbD,permKbD] = cholesky_factorization(KK{nmesh});                      Uinfo.t_chol=toc;
                                                                                tmp=whos('KbD');Uinfo.MB_KbD=tmp.bytes/1048576;
        fprintf(' Factorization of coarse stiffness matrix    = %4.1f sec (%1iMB)\n',Uinfo.t_chol,round(Uinfo.MB_KbD));
    case 'agmg'
        KbD = KK{nmesh}; permKbD = [];
    end
end


% =========================================================================
% PREPARATION OF PRESSURE PRECONDITIONER
% =========================================================================
MM_mu = MM_mu + tril(MM_mu,-1)';
if PC_Pat==1 || OPTS_PC.PC==1
    DMM = diag(MM_mu); % diagonal of MM_mu
    if nsd>1; DMM=sumSDB_fh(DMM,NB,myPnodSDB); end
    DiMM = 1./DMM;  % inverse diagonal of MM_mu
end
if ismember(PC_Pat,[2 3 4]) || ismember(OPTS_PC.PC,[2 3 4])
    LMM = sum(MM_mu,2); % lumped MM_mu
    if nsd>1; LMM=sumSDB_fh(LMM,NB,myPnodSDB); end
    LiMM = 1./LMM;  % inverse lumped MM_mu
end

% if numlabs>1
%     a        = MESH.GCOORD(1,:)';
%     dot_test = COMM.dot(a,a,COMM.unique_nodes{1})
%     a        = MESH.GCOORD(:);
%     dot_test = COMM.dot(a,a,nod2dof(COMM.unique_nodes{1},3))
%     a        = a(ifree);
%     dot_test = COMM.dot(a,a,COMM.unique_Udof{1})
% else
%     a        = MESH.GCOORD(1,:)';
%     dot_test = dot(a,a)
%     a        = MESH.GCOORD(:);
%     dot_test = dot(a,a)
%     a        = a(ifree);
%     dot_test = dot(a,a)
% end

% if numlabs>1
%     save([COMM.prefix '_Workspace3']);
% else
%     WS_Lab02x01 = load('Lab02x01_Workspace3');
%     WS_Lab02x02 = load('Lab02x02_Workspace3');
% end


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
if rmPnullspace; r = rm0space_p(r); end
rP_rms  = normdf_p(r); % norm of residual vector (shown in convergence plot)
tol_Pat = rtol_Pat*rP_rms;
tol_U   = rtol_U*tol_Pat;


% % Show profiling data in terminal
% % ===============================
% fprintf(' Absolute tolerances to calculate initial rP = %0.2e\n',tol_U0);
% fprintf(' Relative/absolute tolerances for pressure   = %0.2e / %0.2e\n',rtol_Pat,tol_Pat);
% fprintf(' Relative/absolute tolerances for velocity   = %0.2e / %0.2e\n',rtol_U,tol_U);

ind = find(abs(U0)>1e-12*max(abs(U0)));
if isempty(ind)
    dU = max(abs(U));
else
    dU = max(abs(U(ind)./U0(ind)));
end
if nsd>1; dU=COMM.maxLabs(dU); end
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
        beta = dot_p(r-rlast,d)/rdlast;  % Polak-Ribiere version 1
%         beta = dot_p(d-dlast,r)/rdlast;  % Polak-Ribiere version 2
        q    = d + beta*q; % calculate NEW search direction
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
    y  = GG*q;                                                              time_Pat=time_Pat+toc;tic
    [z(ifree),nit,time,rU_rms_vec] = mgpcg_solver_p...
            (K1,KK,KbD,permKbD,y(ifree),z(ifree),Ic2f_U,nmg,COMM,OPT_mgpcg);   time_U=time_U+time;nitU=nitU+nit;tic
    Sq = GG'*z;
    if nsd>1; Sq = sumSDB_fh(Sq,NB,myPnodSDB); end
    
    qSq   = dot_p(q,Sq); % denominator in calculating alpha
    rlast = r;      % needed for Polak-Ribiere version 1
%     dlast = d;      % needed for Polak-Ribiere version 2
    alpha = rd/qSq; % step size in direction q
    
    % Update solution and residual
    % ============================
    P = P + alpha*q;  % update pressure solution
    U = U + alpha*z;  % update velocity by accumulating alpha*z
    r = r - alpha*Sq; % update residual vector
    if rmPnullspace; r = rm0space_p(r); end % remove nullspace from both q AND Sq

    rP_rms = normdf_p(r); % norm of residual vector (shown in convergence plot)
    if FigNo
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
    if nsd>1; dU=COMM.maxLabs(dU); end
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
            if rmPnullspace; r = rm0space_p(r); end
            rP_rms   = normdf_p(r); % norm of residual vector (shown in convergence plot)
            
            restart_CG  = 1;
            tol_Pat  = rtol_Pat*rP_rms;
            tol_U    = rtol_U*tol_Pat;
        end
    end
end % END OF PRESSURE ITERATIONS

% if rP_rms>tol_Pat
%     fprintf('\n WARNING: Pressure iterations did not converge within %3i iterations: tol=%0.2e rP_rms=%0.2e\n',...
%             itPat,tol_Pat,rP_rms);
% end
fprintf('\n Time spent on pressure problem               = %6.1f sec\n',time_Pat);
fprintf(' Time spent on velocity problem               = %6.1f sec\n',time_U);
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
fprintf(' Time spent on pressure preconditioning       = %6.1f sec\n',time_z);
fprintf(' Time spent on S-inverse calculation          = %6.1f sec\n',time_PCPat);
end
fprintf(' Total time for solving flow problem          = %6.1f sec\n',etime(clock,clock0));
fprintf('\n');

if FigNo
    saveas(gcf,[SETTINGS.outdir '/' COMM.prefix '_StokesFlowSolverConvergence'],'png');
end

% Check solution
if any(isnan(U)) || any(isnan(P))
	error('Solution containes NaNs !!!!');
end

VAR.Ur = Fb(1:3:end); % radial component
VAR.Ue = Fb(2:3:end); % East component
VAR.Un = Fb(3:3:end); % North component
U      = RR'*U; % rotate velocity solution back to xyz-coordinates
VAR.Ux = U(1:3:end); % x-velocity solution
VAR.Uy = U(2:3:end); % y-velocity solution
VAR.Uz = U(3:3:end); % z-velocity solution
VAR.P  = P;          % pressure solution

% =========================================================================
% ADJUST MESH FOR SURFACE MOTION
% =========================================================================
if strcmp(SETTINGS.top_surface,'free')
    FigNo      = 0;
    [MESH,VAR] = adjust_mesh_fs(VAR,MESH,COMM,dt,SETTINGS,PHYSICS,NUMSCALE,FigNo);
end

% compare_serial_parallel_3d(MESH,COMM,U,'U',0);
% compare_serial_parallel_3d(MESH,COMM,P,'P',1);

% Save profiling data
% *******************
Uinfo.tP     = time_Pat;
Uinfo.tU     = time_U;
byte2MB      = 1/1048576;
KK1          = KK{1}; %#ok<NASGU>
tmp          = whos('KK1'); clear KK1
Uinfo.memKK  = tmp.bytes * byte2MB;
LL           = KK{nmg}; %#ok<NASGU>
tmp          = whos('LL'); clear LL
Uinfo.memLL  = tmp.bytes * byte2MB;
tmp          = whos('GG');
Uinfo.memGG  = tmp.bytes * byte2MB;
Uinfo.nel    = MESH.nel;
Uinfo.nP     = MESH.nVnod;
Uinfo.nU     = nU;
Uinfo.nfree  = length(ifree);
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
                           'is_absolute_tolerance',0,...
                           'coarse_solver',coarse_solver);
    else
        OPT_mgpcg = struct('useMG',useMG,...
                           'itmax',itmax_U,...
                           'tol'  ,tol_U,...
                           'is_absolute_tolerance',1,...
                           'coarse_solver',coarse_solver);
    end
    
    [U(ifree),nit0,tt0,rms_vec] = mgpcg_solver_p...
        (K1,KK,KbD,permKbD,Rhs,U(ifree),Ic2f_U,nmg,COMM,OPT_mgpcg);         time=time+tt0;nit=nit+nit0;rU_rms_vec(end+1:end+length(rms_vec))=rms_vec;tic

    % r = -rhs_P - S p + Fpdil
    %   = -G'*Kinv*Fb - G'*Kinv*G*p + Fpdil , but Fp=G*p and U=Kinv*F (see above)
    % 	= -G'*u_ib - G'*Kinv*Fp + Fpdil     , F=Fb+Fp
    % 	= -G'*U + Fpdil
    r     = -GG'*U + Fpdil; % initial pressure residual vector
    if nsd>1; r = sumSDB_fh(r,NB,myPnodSDB); end
    if rmPnullspace; r = rm0space_p(r); end
    rP_rms  = normdf_p(r); % norm of residual vector (shown in convergence plot)
    tol_Pat = rtol_Pat*rP_rms;
    tol_U   = rtol_U*tol_Pat;

    
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
    
    if tol_U>rU_rms_vec(end)
%         fprintf(' Initial pressure residual was calculated.\n');
        break
    end
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
                if nsd>1; Md = sumSDB_fh(Md,NB,myPnodSDB); end
                d  = d + LiMM .* (r - Md);
            end
        case 4 % inexact Patera algorithm solving S d = r
            OPTS_PC.tol = OPTS_PC.rtol * rP_rms;
            OPTS_PC.tol = max(OPTS_PC.tol,0.75*tol_Pat);
            [d,z,profiling] = inexact_Patera_p(r,OPTS_PC);
    end
end % END OF NESTED FUNCTION precondition_P

% =========================================================================

function display_profiling_data(flag)
if flag==0
    fprintf(' Profiling of iterative velocity-pressure solution...\n\n');
    if PC_Pat==4
        fprintf('                     Patera algorithm                   ||      inexact Patera algorithm     \n');
        fprintf(' -------------------------------------------------------||---------------------------------- \n');
        fprintf('  itP |  time   |    rP    |    rU    |    dU    | nitU ||    rP    | nitP |    rU    | nitU  \n');
        fprintf(' %4i | %7.1f | %0.2e | %0.2e | %0.2e | %4i ||          |      |          |       \n',...
                0,etime(clock,clock0),rP_rms,rU_rms_vec(end),dU,nit);
    else
        fprintf(' Profiling of iterative velocity-pressure solution...\n\n');
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

function [d2,z0,profiling] = inexact_Patera_p(rr,OPTS_PC)

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
    
    rr_rms        = normdf_p(rr);
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
            beta2 = dot_p(rr-rlast2,dd)/rd2last;  % Polak-Ribiere version 1
%             beta2 = dot_p(dd-dlast2,rr)/rd2last; % Polak-Ribiere version 2
            q2    = dd + beta2*q2; % calculate NEW search direction
        end

        if rmPnullspace; q2 = rm0space_p(q2); end % *SUBFUNCTION*
        rd2 = dot_p(rr,dd); % numerator for calculating alpha

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
        [z(ifree),itz1,ttz1,rz_rms_vec] = mgpcg_solver_p...
                (K1,KK,KbD,permKbD,y(ifree),z(ifree),Ic2f_U,nmg,COMM,OPT2_mgpcg); time_z=time_z+ttz1;itz=itz+itz1;tic
        Sq = GG'*z;
        if nsd>1; Sq = sumSDB_fh(Sq,NB,myPnodSDB); end
        
        qSq    = dot_p(q2,Sq); % denominator in calculating alpha
        rlast2 = rr;      % needed for Polak-Ribiere version 1
%         dlast2 = d2;      % needed for Polak-Ribiere version 2
        alpha  = rd2/qSq; % step size in direction q

        % Update solution and residual
        % ============================
        d2 = d2 + alpha*q2; % update pressure "defect"
        rr = rr - alpha*Sq; % update residual of "S d = r"
        z0 = z0 + alpha*z;  % update velocity correction by accumulating alpha*z
        if rmPnullspace; rr = rm0space_p(rr); end % remove nullspace from both q AND Sq

        rr_rms = normdf_p(rr); % norm of residual vector (shown in convergence plot)
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
%         fprintf(' Total time for solving flow problem = %4.1f sec\n',etime(clock,clock0)); %...
%         fprintf(' Total number of Patera iterations   = %4i\n',it_PCPat);
%         fprintf(' Total number of velocity iterations = %4i\n',itz);
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
                    if nsd>1; Md = sumSDB_fh(Md,NB,myPnodSDB); end
                    dd = dd + LiMM .* (rr - Md);
                end
        end
    end % END OF NESTED FUNCTION precondition_d
    
end % END OF NESTED FUNCTION inexact_Patera_p

% =========================================================================

function rms = normdf_p(a)
    if nsd>1
        rms = dot_fh(a,a,unique_Pdofs);
        rms = sqrt(rms/nP_D);
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

end % END OF FUNCTION flow3d_patera_mg_inexact_sph_p

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function [KbD,permKbD,Uinfo] = merge_KbD_v1(Kb,COMM,coarse_solver)

myid = COMM.myid;
nsd  = COMM.nsd;

% Get sparse pattern of BASE stiffness matrix Kb in this subdomain
% "Kbi" and "Kbj" are the locations of the nonzero values "Kbv"
% within the sparse matrix Kb of coarsest mesh (BASE)
[Kbi,Kbj,Kbv] = find(tril(Kb)); % only lower half of Kc is needed
Kbi           = int32(Kbi);
Kbj           = int32(Kbj);
nSDdofs       = length(Kbi);    % number of non-zero entries

% Sum up the number of non-zeros of all subdomains to allocate memory
nDdofs   = COMM.sum_all(nSDdofs);

% Allocate memory for triple format indices of whole-domain coarse KK
KbiD_all = zeros(nDdofs,1,'int32');
KbjD_all = zeros(nDdofs,1,'int32');
KbvD_all = zeros(nDdofs,1);

% Now broadcast the triplet of the subdomain Kb (the pointer
% myUdof_SD2D is used to translate the subdomain dofs into domain dofs)
myUdof_SD2D = COMM.Udof_SD2D_c{myid};
ii          = 1;                                                            MB_Kb=0;tic
for isd=1:nsd
    labBarrier
    if myid==isd
        % I broadcast my domain dofs (rows)
        KbiD = labBroadcast( myid , myUdof_SD2D(Kbi) );
    else
        % I receive neighbor domain dofs (rows)
        KbiD = labBroadcast( isd );
    end
    labBarrier
    if myid==isd
        % I broadcast my domain dofs (columns)
        KbjD = labBroadcast( myid , myUdof_SD2D(Kbj) );
    else
        % I receive neighbor domain dofs (columns)
        KbjD = labBroadcast( isd );
    end
    labBarrier
    if myid==isd
        % I broadcast my coarse stiffness matrix
        KbvD = labBroadcast( myid , Kbv ); clear KKtmp
    else
        % I receive neighbor coarse stiffness matrices
        KbvD = labBroadcast( isd );
    end

    % Construct the triplet of the whole domain coarse stiffness matrix
    if all(KbiD==KbjD)
        error('KbiD and KbjD are identical! Communication screwed up in "merge_Kc"?')
    end
    n                   = length(KbiD); % same as length(jDdofs)
    KbiD_all(ii:ii+n-1) = KbiD(:);
    KbjD_all(ii:ii+n-1) = KbjD(:);
    KbvD_all(ii:ii+n-1) = KbvD(:);
    ii                  = ii + n;
                                                                            tmp=whos('KbvD');MB_Kb=MB_Kb+tmp.bytes/1048576;
                                                                            tmp=whos('KbiD');MB_Kb=MB_Kb+2*tmp.bytes/1048576;
end
% Get rid of arrays that are no longer required
clear Kbi Kbj Kbv KbiD KbiD KbjD KbvD
                                                                            Uinfo.t_broadcast_Kb=toc;Uinfo.MB_Kb=MB_Kb;
fprintf(' Broadcasting coarse stiffness matrices      = %4.1f sec (%1iMB)\n',Uinfo.t_broadcast_Kb,round(MB_Kb));

% Create the DOMAIN BASE stiffness matrix KbD using the above triplet
KbD = sparse3(KbiD_all,KbjD_all,KbvD_all);                                  tmp=whos('KbD');Uinfo.MB_KbD=tmp.bytes/1048576;Uinfo.t_KbD=toc;
fprintf(' Assembly of global coarse stiffness matrix  = %4.1f sec (%1iMB)\n',Uinfo.t_KbD,round(Uinfo.MB_KbD));

% Get rid of arrays that are no longer required
clear KbiD_all KbjD_all KbvD_all 

if strcmp(coarse_solver,'chol')
% Perform Cholesky factorization of DOMAIN stiffness matrix on coarses mesh
% =========================================================================
[KbD,permKbD] = cholesky_factorization(KbD);                                Uinfo.t_chol=toc;tmp=whos('KbD');Uinfo.MB_KbD=tmp.bytes/1048576;
fprintf(' Factorization of coarse stiffness matrix    = %4.1f sec (%1iMB)\n',Uinfo.t_chol,round(Uinfo.MB_KbD));
else
permKbD = [];
end

end % END OF SUBFUNCTION merge_KbD_v1

% #########################################################################

function [KbD,permKbD,Uinfo] = merge_KbD_v2(Kb,COMM,coarse_solver)

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
    fprintf(' Receiving coarse subdomain stiffness matrices = %4.1f sec (%1iMB)\n',Uinfo.t_broadcast_Kb,round(MB_Kb));

    % Create the DOMAIN BASE stiffness matrix KbD using the above triplet
    KbD = sparse3(KbiD_all,KbjD_all,KbvD_all);                              tmp=whos('KbD');Uinfo.MB_KbD=tmp.bytes/1048576;Uinfo.t_KbD=toc;
    % Get rid of arrays that are no longer required
    clear KbiD_all KbjD_all KbvD_all
    fprintf(' Assembly of global coarse stiffness matrix    = %4.1f sec (%1iMB)\n',Uinfo.t_KbD,round(Uinfo.MB_KbD)); 

    if strcmp(coarse_solver,'chol')
    % Perform Cholesky factorization of DOMAIN stiffness matrix on coarses mesh
    % =========================================================================
    [KbD,permKbD] = cholesky_factorization(KbD);                            Uinfo.t_chol=toc;tmp=whos('KbD');Uinfo.MB_KbD=tmp.bytes/1048576;
    fprintf(' Factorization of coarse stiffness matrix      = %4.1f sec (%1iMB)\n',Uinfo.t_chol,round(Uinfo.MB_KbD));
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