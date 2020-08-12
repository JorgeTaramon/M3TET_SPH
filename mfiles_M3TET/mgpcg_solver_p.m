function [x,itCG,time,rrms_vec] = mgpcg_solver_p...
             (K1,K,KbD,permKbD,b,x,Ic2f_U,nmg,COMM,OPTS)
% Usage: [x,itCG,time,rrms_vec] = mgpcg_solver_p...
%            (K1,K,LbD,permLbD,b,x,Ic2f_U,nmg,COMM,OPTS)
%
% Purpose: Solution algorithm for the equation K{1} x = b
%          Method: Conjugate Gradient algorithm that is preconditioned by a
%          single V-cycle of geometric multigrid; direct Cholesky solver on
%          coarsest multigrid level.
% Input:
%   K1      :: coefficient matrix for finest MG level
%   K       :: cell array storing coefficient matrices all but the finest
%              MG level (K{1} is used for diagonal of K1)
%   KbD     :: Cholesky factorization of coarsest level K{nmesh}
%   permKbD :: premutation vector for Cholesky solver on coarsest MG level
%   b       :: right-hand-side vector
%   x       :: guess for the solution vector
%   Ic2f_U  :: cell array storing interplation matrices to transfer data
%              between neighboring MG levels (coarser to finer)
%   nmg     :: number of geometric multigrid (MG) levels
%   COMM    :: structure containing all SD communication data
%   OPTS    :: structure containing control parameters for the iterative
%              solver
% Output:
%   x        :: solution vector
%   itCG     :: number of CG iterations until convergence
%   time     :: computer time
%   rrms_vec :: norm of residual vector in each CG iteration
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Feb 2011
% JH May 2014: - included sparse-matrix-vector-multiplication from MUTILS
%              - added Yvan Notay's algebraic multigrid solver AGMG as
%                alternative solver on the coarsest mesh (it is slower than
%                Cholesky if the mesh is too small)
% JH Mar 2015: cleaned up
%
t0   = tic;

myid = COMM.myid;
nsd  = COMM.nsd;
if nsd>1
    NB        = COMM.NB;      % list of neighboring SDs
    sumSDB_fh = COMM.sum_SDB; % function handle
    nUD       = COMM.nUD;     % number of free vel dof in whole domain on
                              % each multigrid level
    myUdofSDB   = COMM.myUdofSDB; % subdomain boundary information for all
                                  % neighbors on all multigrid levels
end

check_input_args(); % *NESTED FUNCTION*
nmesh = length(K);  % number of meshes (is equal to nmg if no extra-coarse 
                    % BASE mesh is used as the lowest MG level).

[useMG,rtol,tol,itmax,FigNo,norm_res_vs_time,PC,Jwght_down,Jwght_up,...
    coarse_solver] = setup_mgpcg(OPTS,nmg,nmesh); % *SUBFUNCTION*

if norm_res_vs_time
    time = zeros(itmax+1,1);
end

% Extract diagonal from coefficient matrices, invert and store in cell
% array (needed for Jacobi smoother)
KDi = cell(nmesh,1);
for jmg=1:nmesh
    if jmg==1
        KD = K{1};
    else
        KD = full(diag(K{jmg})); % diagonal of K{jmg}
    end
    if nsd>1
        KD = sumSDB_fh(KD,NB,myUdofSDB{jmg});
    end
    KDi{jmg} = 1./KD;        % inverse diagonal of K{jmg}
end


% Calculate initial residual
% ==========================
r  = b - mat_x_vec(K1,x); % *SUBFUNCTION*
if nsd>1
    r = sumSDB_fh(r,NB,myUdofSDB{1});
end
rrms = normdf_p(r,1); % *NESTED FUNCTION

if tol==0
    if rtol<=0
        error('Relative tolerance must be a positive value!!');
    end
    tol = rtol*rrms;  % tolerance for convergence
end

rrms_vec    = zeros(itmax+1,1);
rrms_vec(1) = rrms;

if FigNo
    figure(FigNo);clf
    if norm_res_vs_time
        plot(toc(t0),log10(rrms),'k.'); set(gca,'FontSize',14);
        xlabel('Time (sec)'); ylabel('Norm of residual vector');
    else
        plot(0,log10(rrms),'k.'); set(gca,'FontSize',14);
        xlabel('Iteration'); ylabel('Norm of residual vector');
    end
    hold on
    drawnow
end

% Begin of CG iterations
% ======================
for it=1:itmax
    
    % Preconditioning (approximate error of current x)
    % ================================================
    if PC==3 % multigrid preconditioner (best)
%         d = MG_V_cycle_v1;   % *NESTED FUNCTION*
        d = MG_V_cycle_v2;   % *NESTED FUNCTION*
    else     % other preconditioners
        d = precondition; % *NESTED FUNCTION*
    end
    
    if useMG
        x              = x + d;
        r              = b - mat_x_vec(K1,x); % *SUBFUNCTION*
        if nsd>1; r = sumSDB_fh(r,NB,myUdofSDB{1}); end
        rrms           = normdf_p(r,1);
        rrms_vec(it+1) = rrms;
        if FigNo
            figure(FigNo);
            if norm_res_vs_time
                plot(toc(t0),log10(rrms),'k.');
            else
                plot(it,log10(rrms),'k.');
            end
            drawnow
        end
        if norm_res_vs_time; time(it+1) = toc(t0); end
        if rrms<tol
            break
        else
            continue
        end
    end
    
    % get new search direction q
    % ==========================
    if it==1
        % Initialize Conjugate Gradients
        % ==============================
        q  = d; % define FIRST search direction q

    else
        rdlast = rd; % save old rd (for calculating beta)
        % Make new search direction q K-orthogonal to all previous q's
        beta = dot_p(r-rlast,d,1)/rdlast;  % Polak-Ribiere version 1
%         beta = dot_p(d-dlast,r,1)/rdlast;  % Polak-Ribiere version 2
        q    = d + beta*q; % calculate NEW search direction
    end

    rd    = dot_p(r,d,1);    % numerator for calculating alpha
    Kq    = mat_x_vec(K1,q); % calculate K1*q, *SUBFUNCTION*
    if nsd>1; Kq=sumSDB_fh(Kq,NB,myUdofSDB{1}); end
    qKq   = dot_p(q,Kq,1);   % denominator in calculating alpha
    
    rlast = r;         % needed for Polak-Ribiere version 1
%     dlast = d;         % needed for Polak-Ribiere version 2
    alpha = rd/qKq;    % step size in direction q

    % Update solution and residual
    % ============================
    x = x + alpha * q;  % update solution
    r = r - alpha * Kq; % update residual
    
    rrms           = normdf_p(r,1); % norm of residual vector
    rrms_vec(it+1) = rrms;
    if norm_res_vs_time; time(it+1) = toc(t0); end
    
    if FigNo
        figure(FigNo);
        if norm_res_vs_time
            plot(toc(t0),log10(rrms),'k.');
        else
            plot(it,log10(rrms),'k.');
        end
        drawnow
    end
    
    % Check if solution converged
    % ===========================
    if rrms<tol
        break
    end
end % end of CG loop

if norm_res_vs_time
    time(it+2:end) = [];
else
    time = toc(t0);
end
rrms_vec(it+2:end) = [];
itCG = it;

if FigNo
    figure(FigNo);clf
    if norm_res_vs_time
        plot(time,log10(rrms_vec),'k.-'); set(gca,'FontSize',14);
        xlabel('Time (sec)'); ylabel('Norm of residual vector');
    else
        plot(log10(rrms_vec),'k.-'); set(gca,'FontSize',14);
        xlabel('Iteration'); ylabel('Norm of residual vector');
    end
end

% =========================================================================
%                            NESTED FUNCTIONS
% =========================================================================

function rms = normdf_p(a,jmg)
    if nsd>1
        rms = COMM.normdf(a,COMM.unique_Udof{jmg});
    else
        rms = normdf(a);
    end
end % END OF NESTED FUNCTION normdf_p

% =========================================================================

function ab = dot_p(a,b,jmg)
    if nsd>1
        ab = COMM.dot(a,b,COMM.unique_Udof{jmg});
    else
        ab = dot(a,b);
    end
end % END OF NESTED FUNCTION dot_p

% =========================================================================

function d = precondition
    switch PC
        case 0 % No preconditioning (i.e. CG algorithm)
            d = r;
            
        case  1 % Diagonal scaling (Jacobi)
            d = r.*KDi{1};

        case 2 % Symmetric Gauss-Seidel
            if nsd>1; error('Gauss-Seidel preconditioner does not work in parallel!'); end
            w = tril(K1) \ r;  % FORWARD SOLVE
            w = w .* K{1};     % multiply diagonal
            d = tril(K1)' \ w; % BACKWARD SOLVE
    end
end % END OF NESTED FUNCTION precondition

% =========================================================================

function d = MG_V_cycle_v1
    % Performs a single multigrid V-cycle.
    % Returns defect vector d for changing input r on finest MG level.
    r0    = r;             % keep original residual r
    dd    = cell(nmesh-1,1);
    rr    = cell(nmesh-1,1);
    
    % Pre-smoothing of residual (downward path of V-cycle)
    % ====================================================
    for img=1:nmesh-1 % downward loop over MG-level

        rr{img} = r0; % Note: r0 changes in size and value after each MG
                      %       level as it is restricted
        r1      = r0;
        nU_iMG  = length(r0);
        d1      = zeros(nU_iMG,1);
        % Relaxations (smoothing) on MG level "img"
        for ii = 1:length(Jwght_down{img})
            d1 = d1 + Jwght_down{img}(ii)*r1.*KDi{img};  % Jacobi (diagonal) smoothing
            if img==1
                Kd = mat_x_vec(K1,d1); % *SUBFUNCTION*
            else
                Kd = mat_x_vec(K{img},d1); % *SUBFUNCTION*
            end
            if nsd>1; Kd=sumSDB_fh(Kd,NB,myUdofSDB{img}); end % *NESTED FUNCTION*
            r1 = r0 - Kd; % calculate new residual
        end
        
        dd{img} = d1; % keep defect vector for this level
                      % (used later during the upward path)
        
        if nsd>1
            % Parallel: r1 is summed up along SD boundaries (because r0 and Kd
            % are summed up, so is r1 = r0 - Kd)
            % Before restricting r1, the non-unique dofs must be set to zero to
            % avoid multiple contributions of dofs at SD boundaries.
            % nuniqUdofs = setdiff(1:length(r1),uniqUdofs{img});
            i0     = setdiff(1:nU_iMG,COMM.unique_Udof{img});
            r1(i0) = 0;
        end
        
        % Restrict r1 to the next coarser level
        r0 = Ic2f_U{img}' * r1;
        
        % Sum up residual along SD boundary (not if we arrived on the
        % coarsest MG level)
        if nsd>1 && img<nmesh-1
            r0=sumSDB_fh(r0,NB,myUdofSDB{img+1}); % *NESTED FUNCTION*
        end
    end

    % Solve K{nmesh} d = r on the coarsest MG level
    % =============================================
    if nsd>1
        % Broadcast SD residuals to form domain residual
        r_c = accumulate_domain_residual_v1();
    else
        r_c = r0;
    end
    
    % Note: KbD is the Cholesky factorization of the merged K matrices of
    % all subdomains on their coarsest MG level.
    d = cholesky_solve(r_c,KbD,permKbD);
    
    % Extract the subdomain part of the domain solution "d"
    if nsd>1
        d = d( COMM.Udof_SD2D_c{myid} );
    end

    
    % Post-smoothing of residual (upward path of V-cycle)
    % ===================================================
    for img=nmesh-1:-1:1 % upward loop over MG-levels
        % Interpolate coarse error to next finer level and add to error
        % component "dd" that was saved after each downward level
        di = Ic2f_U{img}*d(:);
        
        % The next summation along SD boundaries is only necessary after
        % interpolating the error from the BASE mesh to the coarsest
        % geometric MG mesh
        if nsd>1 && nmesh>nmg && img==nmesh-1
            di=sumSDB_fh(di,NB,myUdofSDB{img}); % *NESTED FUNCTION*
        end

        d = dd{img} + di;
        
        % Relaxations (smoothing) on MG level "img"
        for ii = 1:length(Jwght_up{img})
            if img==1
                Kd = mat_x_vec(K1,d); % *SUBFUNCTION*
            else
                Kd = mat_x_vec(K{img},d); % *SUBFUNCTION*
            end
            if nsd>1; Kd=sumSDB_fh(Kd,NB,myUdofSDB{img}); end % *NESTED FUNCTION*
            r1 = rr{img} - Kd;
            d  = d + Jwght_up{img}(ii)*r1.*KDi{img}; % Jacobi (diagonal) smoothing
        end
    end

    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    %                            NESTED FUNCTIONS
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    
    function r_c = accumulate_domain_residual_v1
        % Broadcast SD residuals to form domain residual
        if myid==1
            r_c = zeros(nUD(nmesh),1); % initialize domain residual
            r_c(COMM.Udof_SD2D_c{1}) = r0;  % write residual of lab 1
            for isd=2:nsd
                % now receive residuals from labs 2:nsd and accumulate
                UdofsD = COMM.Udof_SD2D_c{isd};
                r_c(UdofsD) = r_c(UdofsD) + labReceive( isd );
            end
            labBarrier
            labBroadcast( myid , r_c ); % send the domain residual to all other labs
        else
            labSend( r0 , 1 ); % send SD residual to lab 1
            labBarrier
            r_c = labBroadcast( 1 ); % receive domain residual from lab 1
        end
    end % END OF NESTED FUNCTION accumulate_domain_residual
    
end % END OF NESTED FUNCTION MG_V_cycle_v1

% =========================================================================

function d = MG_V_cycle_v2
    % Performs a single multigrid V-cycle.
    % Returns defect vector d for changing input r on finest MG level.
    r0    = r;             % keep original residual r
    dd    = cell(nmesh-1,1);
    rr    = cell(nmesh-1,1);
    
    % Pre-smoothing of residual (downward path of V-cycle)
    % ====================================================
    for img=1:nmesh-1 % downward loop over MG-level

        rr{img} = r0; % Note: r0 changes in size and value after each MG
                      %       level as it is restricted
        r1      = r0;
        nU_iMG  = length(r0);
        d1      = zeros(nU_iMG,1);
        % Relaxations (smoothing) on MG level "img"
        for ii = 1:length(Jwght_down{img})
            d1 = d1 + Jwght_down{img}(ii)*r1.*KDi{img};  % Jacobi (diagonal) smoothing
            if img==1
                Kd = mat_x_vec(K1,d1); % *SUBFUNCTION*
            else
                Kd = mat_x_vec(K{img},d1); % *SUBFUNCTION*
            end
            if nsd>1; Kd=sumSDB_fh(Kd,NB,myUdofSDB{img}); end % *NESTED FUNCTION*
            r1 = r0 - Kd; % calculate new residual
        end
        
        dd{img} = d1; % keep defect vector for this level
                      % (used later during the upward path)
        
        if nsd>1
            % Parallel: r1 is summed up along SD boundaries (because r0 and Kd
            % are summed up, so is r1 = r0 - Kd)
            % Before restricting r1, the non-unique dofs must be set to zero to
            % avoid multiple contributions of dofs at SD boundaries.
            % nuniqUdofs = setdiff(1:length(r1),uniqUdofs{img});
            i0     = setdiff(1:nU_iMG,COMM.unique_Udof{img});
            r1(i0) = 0;
        end
        
        % Restrict r1 to the next coarser level
        r0 = Ic2f_U{img}' * r1;
        
        % Sum up residual along SD boundary (not if we arrived on the
        % coarsest MG level)
        if nsd>1 && img<nmesh-1
            r0=sumSDB_fh(r0,NB,myUdofSDB{img+1}); % *NESTED FUNCTION*
        end
    end

    
    % Solve K{nmesh} d = r on the coarsest MG level
    % =============================================
    switch coarse_solver
    case 'chol'
        if nsd>1
            % Broadcast SD residuals to form domain residual
            r_c = accumulate_domain_residual_v2(); % *NESTED FUNCTION*
            if myid==1
                d_c = cholesky_solve(r_c,KbD,permKbD);
            end
            % Extract the subdomain part of the domain solution "d"
            d = distribute_subdomain_errors(); % *NESTED FUNCTION*

        else
            d = cholesky_solve(r0,KbD,permKbD);
        end
        
    case 'ichol'
        if nsd>1
            error(' Coarse solver "ichol" not yet parallel.');
            % Broadcast SD residuals to form domain residual
            r_c = accumulate_domain_residual_v2(); % *NESTED FUNCTION*
            if myid==1
                d_c = cholesky_solve(r_c,KbD,permKbD);
            end
            % Extract the subdomain part of the domain solution "d"
            d = distribute_subdomain_errors(); % *NESTED FUNCTION*

        else
            d = cholesky_solve(r0,KbD,permKbD);
        end
        
        % Relaxations (smoothing) on MG level "nmg"
        Jwght_nmg = [0.3 0.6 1 0.6 0.3];
        for ii = 1:length(Jwght_nmg)
            Kd = mat_x_vec(K{nmg},d); % *SUBFUNCTION*
            if nsd>1; Kd=sumSDB_fh(Kd,NB,myUdofSDB{nmg}); end % *NESTED FUNCTION*
            r1 = r0 - Kd;
            d  = d + Jwght_nmg(ii)*r1.*KDi{nmg}; % Jacobi (diagonal) smoothing
        end
            
    case 'agmg'
        icg     = 1;
            % 1, use CG and performs some simplifications based on
            % the assumption that thecoefficient matrix A is symmetric;
            % should be used only when A is symmetric and positive definite.
            % >=2, use GCR restarted each RESTART iterations.
            % 0 or [], then agmg use the default, which is
            % GCR restarted each 10 iterations.
        verbose = 1;
        rtol    = 1e-2;
        itmax   = 50;
        if nsd>1
            % Broadcast SD residuals to form domain residual
            r_c = accumulate_domain_residual_v2(); % *NESTED FUNCTION*
            if myid==1
                [d_c,flag,relres,nit0,rms_vec] = agmg...
                    (K{end},r_c,icg,rtol,itmax,verbose);
            end
            % Extract the subdomain part of the domain solution "d"
            d = distribute_subdomain_errors(); % *NESTED FUNCTION*
        else        
            [d,flag,relres,nit0,rms_vec] = agmg...
                    (K{end},r0,icg,rtol,itmax,verbose);
        end
        
    case 'bslash'
        d = K{end} \ r0;
        
    case 'jacobi'
        % Jacobi iterations (a very coarse mesh is required here !!)
        d   = r0 .* KDi{end};
        Jwght_base = [0.25 0.4 0.9];
        for is=1:length(Jwght_base)
            Kxd = mat_x_vec(K{end},d);
            d   = d + Jwght_base(is) * KDi{end} .* (r0 - Kxd);
        end
        clear Kxd
        
    otherwise
        error(' Unknown coarse solver.');
    end
    
    % Post-smoothing of residual (upward path of V-cycle)
    % ===================================================
    for img=nmesh-1:-1:1 % upward loop over MG-levels
        % Interpolate coarse error to next finer level and add to error
        % component "dd" that was saved after each downward level
        di = Ic2f_U{img}*d(:);
        
        % The next summation along SD boundaries is only necessary after
        % interpolating the error from the BASE mesh to the coarsest
        % geometric MG mesh
        if nsd>1 && nmesh>nmg && img==nmesh-1
            di=sumSDB_fh(di,NB,myUdofSDB{img}); % *NESTED FUNCTION*
        end

        d = dd{img} + di;
        
        % Relaxations (smoothing) on MG level "img"
        for ii = 1:length(Jwght_up{img})
            if img==1
                Kd = mat_x_vec(K1,d); % *SUBFUNCTION*
            else
                Kd = mat_x_vec(K{img},d); % *SUBFUNCTION*
            end
            if nsd>1; Kd=sumSDB_fh(Kd,NB,myUdofSDB{img}); end % *NESTED FUNCTION*
            r1 = rr{img} - Kd;
            d  = d + Jwght_up{img}(ii)*r1.*KDi{img}; % Jacobi (diagonal) smoothing
        end
    end

    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    %                            NESTED FUNCTIONS
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    
    function r_c = accumulate_domain_residual_v2
        % Broadcast SD residuals to form domain residual
        if myid==1
            r_c = zeros(nUD(nmesh),1); % initialize domain residual
            r_c(COMM.Udof_SD2D_c{1}) = r0;  % write residual of lab 1
        else
            r_c = [];
        end
        for isd=2:nsd
            labBarrier
            tag = 1000+isd;
            if myid==1
                % now receive residuals from labs 2:nsd and accumulate
                UdofsD      = COMM.Udof_SD2D_c{isd};
                r_c(UdofsD) = r_c(UdofsD) + COMM.recv_1NB(isd,tag);
            elseif myid==isd
                COMM.send_1NB(1,r0,tag); % send SD residual to lab 1
            end
        end
    end % END OF NESTED FUNCTION accumulate_domain_residual
    
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    function d = distribute_subdomain_errors
        % Extract the subdomain part of the domain solution "d"
        for isd=2:nsd
            labBarrier
            tag = 2000+isd;
            if myid==1
                d_NB = d_c( COMM.Udof_SD2D_c{isd} );
                COMM.send_1NB(isd,d_NB,tag);
            elseif myid==isd
                d = COMM.recv_1NB(1,tag);
            end
        end
        if myid==1
            d = d_c( COMM.Udof_SD2D_c{myid} );
        end
    end % END OF NESTED FUNCTION distribute_subdomain_errors
    
end % END OF NESTED FUNCTION MG_V_cycle_v2

% =========================================================================

function check_input_args()
    % check input
    n = length(b);
    if length(x)~=n
        error('Length of RHS vector does not match solution vector');
    end
    if n~=length(K{1})
        error('Length of RHS vector does not match coefficient matrix');
    end
    if any(isnan(b))
        error('RHS vector continas NaNs!!');
    end
end % END OF NESTED FUNCTION check_input_args

end % END OF FUNCTION mgpcg_solver_p

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function [useMG,rtol,tol,itmax,FigNo,norm_res_vs_time,PC,...
          Jwght_down,Jwght_up,coarse_solver] = setup_mgpcg(OPTS,nmg,nmesh)

% --------------------------- DEFAULT SETTINGS ---------------------------
useMG = 0;    % use multigrid-preconditioned CG
rtol  = 0;    % don't use relative tolerance but...
tol   = 1e-8; % this absolute tolerance
itmax = 50;   % maximum number of iterations
FigNo = 0;            % set =0 to suppress plotting
norm_res_vs_time = 0; % =0 --> plot norm residual vs iteration
                      % =1 --> plot norm residual vs computer time
coarse_solver = 'chol'; % solver on coarsest multigrid level
PC = 3; %  0 ==> no preconditioning (d=r)
        %  1 ==> Jacobi diagonal scaling
        %  2 ==> symmetric Gauss-Seidel
        %  3 ==> single V-cycle of geometric multigrid
        
% Damping factors for Jacobi smoother on every MG level
% (number of smoothing sweeps equal to length of Jwght{jmg})
Jwght_down = cell(nmesh-1,1);
Jwght_up   = cell(nmesh-1,1);

% % * * * * * * * Setup 1 * * * * * * * 
% % Define smoothing weights for FINEST (time intensive) level
% Jwght_down{1} = [0.3 0.6 1];
% Jwght_up{1}   = Jwght_down{1}(end:-1:1); % to keep symmetry
% % Define smoothing weights for ALL OTHER levels
% for jmg=2:nmg-1
%     Jwght_down{jmg} = [0.3 0.6 1];
%     Jwght_up{jmg}   = Jwght_down{jmg}(end:-1:1); % to keep symmetry
% end
% % In case we're using a BASE mesh
% if nmesh~=nmg 
%     Jwght_down{nmesh-1} = [0.3 0.6 1 0.6 0.3];
%     Jwght_up{nmesh-1}   = Jwght_down{nmesh-1}(end:-1:1); % to keep symmetry
% end

% % * * * * * * * Setup 2 * * * * * * * 
% % Define smoothing weights for FINEST (time intensive) level
% Jwght_down{1} = [0.6];
% Jwght_up{1}   = [1 0.6 0.3];
%     % this choice makes it an unsymmetric operation (could cause problems,
%     % BE CAREFUL!!)
% % Define smoothing weights for ALL OTHER levels
% for jmg=2:nmg-1
%     Jwght_down{jmg} = [0.6];
%     Jwght_up{jmg}   = [1 0.6 0.3];
% end
% % In case we're using a BASE mesh
% if nmesh~=nmg 
%     Jwght_down{nmesh-1} = [0.3 0.6 1 0.6 0.3];
%     Jwght_up{nmesh-1}   = Jwght_down{nmesh-1}(end:-1:1); % to keep symmetry
% end

% * * * * * * * Setup 3 * * * * * * * 
% Define smoothing weights for FINEST (time intensive) level
Jwght_down{1} = [0.9 0.4 0.25];
Jwght_up{1}   = [0.25 0.4 0.9];
% Define smoothing weights for ALL OTHER levels
for jmg=2:nmg-1
    Jwght_down{jmg} = [0.25 0.4 0.9];
    Jwght_up{jmg}   = Jwght_down{jmg}(end:-1:1); % to keep symmetry
end
% In case we're using a BASE mesh
if nmesh~=nmg 
    Jwght_down{nmesh-1} = [0.25 0.4 0.9 0.4 0.25];
    Jwght_up{nmesh-1}   = Jwght_down{nmesh-1}(end:-1:1); % to keep symmetry
end
% --------------------------- DEFAULT SETTINGS ---------------------------


% Overwrite default values if provided in struture "OPTS"
if isfield(OPTS,'useMG')
    useMG = OPTS.useMG;
end
if isfield(OPTS,'coarse_solver')
    coarse_solver = OPTS.coarse_solver;
end
if isfield(OPTS,'FigNo')
    FigNo = OPTS.FigNo;
end
if isfield(OPTS,'norm_res_vs_time')
    norm_res_vs_time = OPTS.norm_res_vs_time;
end
if isfield(OPTS,'tol') && isfield(OPTS,'is_abs_tol')
    if OPTS.is_abs_tol
        tol  = OPTS.tol; rtol = 0;
    else
        rtol = OPTS.tol; tol = 0;
    end
end
if isfield(OPTS,'itmax')
    itmax = OPTS.itmax;
end
if isfield(OPTS,'PC')
    PC = OPTS.PC;
end
if isfield(OPTS,'Jwght_down')
    Jwght_down = OPTS.Jwght_down;
end
if isfield(OPTS,'Jwght_up')
    Jwght_up = OPTS.Jwght_up;
end

end % END OF FUNCTION setup_mgpcg

% #########################################################################

function a = mat_x_vec(K,x)

% Times in brackets are for 65 muliplications with nU==225370
if isstruct(K)
    a = spmv(K,x); % (4.1-4.5 sec)
                   % NOTE: no input "opts" if K is in mutils format
else
    if max(abs(K(1,2:end)))==0
        a = K*x + K'*x - full(diag(K)).*x; % performs: a = K * x; (10.1-10.6 sec)
%         a = K*x + tril(K,-1)'*x; % performs: a = K * x; SLOWER (17.5-17.8 sec)
    else
        a = K*x; % (6.4-6.9 sec)
    end
end

end % END OF SUBFUNCTION mat_x_vec