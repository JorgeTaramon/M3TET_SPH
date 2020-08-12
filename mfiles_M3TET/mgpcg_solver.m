function [x,itCG,time,rrms_vec] = mgpcg_solver(K,LbD,permLbD,b,x,Ic2f_U,nmg,OPT)
% Usage: [x,itCG,time,rrms_vec] = mgpcg_solver(K,LbD,permLbD,b,x,Ic2f_U,nmg,OPT)
%
% Purpose: Solution algorithm for the equation K{1} x = b
%          Method: Conjugate Gradient algorithm that is preconditioned by a
%          single V-cycle of geometric multigrid; direct Cholesky solver on
%          coarsest multigrid level.
% Input:
%   K       :: cell array that stores the coefficient matrix for each MG level
%   LbD     :: Cholesky factorization of coarsest level K{nmesh}
%   permLbD :: premutation vector for Cholesky solver on coarsest MG level
%   b       :: right-hand-side vector
%   x       :: guess for the solution vector
%   Ic2f_U  :: cell array storing interplation matrices to transfer data
%              between neighboring MG levels (coarser to finer)
%   nmg     :: number of geometric multigrid (MG) levels
%   OPT     :: structure containing control parameters for the iterative
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

% JH Feb 2011
%

check_input_args(); % *NESTED FUNCTION
nmesh = length(K);  % number of meshes (is equal to nmg if no extra-coarse 
                    % BASE mesh is used as the lowest MG level).

[useMG,rtol,tol,itmax,FigNo,norm_res_vs_time,PC,Jwght_down,Jwght_up] = ...
    setup_mgpcg(OPT,nmg,nmesh); % *SUBFUNCTION*

% Extract diagonal from coefficient matrices, invert and store in cell
% array (needed for Jacobi smoother)
KDi = cell(nmesh-1,1);
for jmg=1:nmesh-1
    KDi{jmg} = 1./full(diag(K{jmg}));
end
                                                                            tic
% Calculate initial residual
% ==========================
r    = b - K{1}*x;
rrms = normdf(r);  % norm of residual vector

if tol==0
    if rtol<=0
        error('Relative tolerance must be a positive value!!');
    end
    tol = rtol*rrms;  % tolerance for convergence
end

rrms_vec    = zeros(itmax+1,1);
rrms_vec(1) = rrms;
if norm_res_vs_time
    time = zeros(itmax+1,1);
end

% Begin of CG iterations
% ======================
for it=1:itmax
    
    % Preconditioning (approximate error of current x)
    % ================================================
    if PC==3 % multigrid preconditioner (best)
        d = mg_V_cycle;   % *NESTED FUNCTION*
    else     % other (simpler) preconditioners
        d = precondition; % *NESTED FUNCTION*
    end
    
    if useMG
        x              = x + d;
        r              = b - K{1}*x;
        rrms           = normdf(r);
        rrms_vec(it+1) = rrms;
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
        beta = dot(r-rlast,d)/rdlast;  % Polak-Ribiere version 1
%         beta = dot(d-dlast,r)/rdlast;  % Polak-Ribiere version 2
        q    = d + beta*q; % calculate NEW search direction
    end

    rd    = dot(r,d);  % numerator for calculating alpha
    Kq    = K{1}*q;    % calculate K1*q
    qKq   = dot(q,Kq); % denominator in calculating alpha
    
    rlast = r;         % needed for Polak-Ribiere version 1
%     dlast = d;         % needed for Polak-Ribiere version 2
    alpha = rd/qKq;    % step size in direction q

    % Update solution and residual
    % ============================
    x = x + alpha * q;  % update solution
    r = r - alpha * Kq; % update residual
    
    rrms           = normdf(r); % norm of residual vector
    rrms_vec(it+1) = rrms;
    if norm_res_vs_time; time(it+1) = toc; end
    
    % Check if solution converged
    % ===========================
    if rrms<tol
        break
    end
end % end of CG loop

if norm_res_vs_time
    time(it+2:end) = [];
else
    time = toc;
end
rrms_vec(it+2:end) = [];
itCG = it;

if FigNo
    figure(FigNo);clf
    if norm_res_vs_time
        semilogy(time,rrms_vec,'k.-'); set(gca,'FontSize',14);
        xlabel('Time (sec)'); ylabel('Norm of residual vector');
    else
        semilogy(rrms_vec,'k.-'); set(gca,'FontSize',14);
        xlabel('Iteration'); ylabel('Norm of residual vector');
    end
end

%%%%%%%%%%%%%%%% NESTED FUNCTIONS %%%%%%%%%%%%%%%%%

function d = precondition
    switch PC
        case 0 % No preconditioning (i.e. CG algorithm)
            d = r;
            
        case  1 % Diagonal scaling (Jacobi)
            d = r.*KDi{1};

        case 2 % Symmetric Gauss-Seidel
            w = tril(K{1}) \ r;  % FORWARD SOLVE
            w = w .* diag(K{1}); % multiply diagonal
            d = tril(K{1})' \ w; % BACKWARD SOLVE
    end
end % END OF NESTED FUNCTION precondition

function d = mg_V_cycle
    % Performs a single multigrid V-cycle.
    % Returns defect vector d for changing input r on finest MG level.
    r0    = r;             % keep original residual r
    dd    = cell(nmesh-1,1);
    rr    = cell(nmesh-1,1);
    
    % Pre-smoothing of residual (downward path of V-cycle)
    % ====================================================
    for img=1:nmesh-1 % downward loop over MG-levels
        rr{img} = r0; % Note: r0 changes in size and value after each MG
                      %       level as it is restricted
        r1      = r0;
        d1      = zeros(size(r0));
        % Relaxations (smoothing) on MG level "img"
        for ii = 1:length(Jwght_down{img})
            d1 = d1 + Jwght_down{img}(ii)*r1.*KDi{img}; % Jacobi (diagonal) smoothing
            r1 = r0 - K{img}*d1; % calculate new residual
        end
        dd{img} = d1; % keep defect vector for this level
                      % (used later during the upward path)

        % Restrict r1 to the next coarser level
        r0 = Ic2f_U{img}' * r1;
    end

    
    
    % Note: LbD is the Cholesky factorization of the merged K matrices of
    % all subdomains on their coarsest MG level.
    % Cholesky: (LbD LbD') d = r
    % step (1) LbD  y  = r   (FORWARD)
    % step (2) LbD' d  = tmp (BACKWARD)
%     d          = zeros(size(r0));
%     y          = LbD  \ r0(permLbD); % FORWARD SOLVE of triangular system
%     d(permLbD) = LbDt \ y;           % BACKWARD SOLVE of triangular system
    d          = cholesky_substitution(r0,LbD,permLbD);
    
    
    % Post-smoothing of residual (upward path of V-cycle)
    % ===================================================
    for img=nmesh-1:-1:1 % upward loop over MG-levels

        % Interpolate coarse error to next finer level and add to error
        % component "dd" that was saved after each downward level
        d = dd{img} + Ic2f_U{img}*d(:);
        % Relaxations (smoothing) on MG level "img"
        for ii = 1:length(Jwght_up{img})
            r1 = rr{img} - K{img}*d;
            d  = d + Jwght_up{img}(ii)*r1.*KDi{img}; % Jacobi (diagonal) smoothing
        end
    end
end % END OF NESTED FUNCTION mg_V_cycle

function check_input_args()
    % check input
    n = length(b);
    if length(x)~=n
        error('Length of RHS vector does not match solution vector');
    end
    if n~=size(K{1},1) || n~=size(K{1},2)
        error('Length of RHS vector does not match coefficient matrix');
    end
    if any(isnan(b))
        error('RHS vector continas NaNs!!');
    end
end % END OF NESTED FUNCTION check_input_args

end % END OF FUNCTION mgpcg_solver

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function [useMG,rtol,tol,itmax,FigNo,norm_res_vs_time,PC,Jwght_down,Jwght_up] = ...
    setup_mgpcg(OPT,nmg,nmesh)

% --------------------------- DEFAULT SETTINGS ---------------------------
useMG = 0;    % use multigrid-preconditioned CG
rtol  = 0;    % don't use relative tolerance but...
tol   = 1e-8; % this absolute tolerance
itmax = 50;   % maximum number of iterations
FigNo = 0;            % set =0 to suppress plotting
norm_res_vs_time = 0; % =0 --> plot norm residual vs iteration
                      % =1 --> plot norm residual vs computer time
PC = 3; %  0 ==> no preconditioning (d=r)
        %  1 ==> Jacobi diagonal scaling
        %  2 ==> symmetric Gauss-Seidel
        %  3 ==> single V-cycle of geometric multigrid
        
% Damping factors for Jacobi smoother on every MG level
% (number of smoothing sweeps equal to length of Jwght{jmg})
Jwght_down = cell(nmesh-1,1);
Jwght_up   = cell(nmesh-1,1);
% Define smoothing weights for FINEST (time intensive) level
Jwght_down{1} = [0.3 0.6 1];
Jwght_up{1}   = Jwght_down{1}(end:-1:1);
% Define smoothing weights for ALL OTHER levels
for jmg=2:nmg-1
    Jwght_down{jmg} = [0.3 0.6 1 0.6 0.3];
    Jwght_up{jmg}   = Jwght_down{jmg}(end:-1:1);
end
% In case we're using a BASE mesh
if nmesh~=nmg 
    Jwght_down{nmesh-1} = [0.3 0.6 1 1.3 1 0.6 0.3];
    Jwght_up{nmesh-1}   = Jwght_down{nmesh-1}(end:-1:1);
end
% --------------------------- DEFAULT SETTINGS ---------------------------

% Overwrite default values if provided in struture "OPT"
if isfield(OPT,'useMG')
    useMG = OPT.useMG;
end
if isfield(OPT,'norm_res_vs_time')
    norm_res_vs_time = OPT.norm_res_vs_time;
end
if isfield(OPT,'tol') && isfield(OPT,'is_absolute_tolerance')
    if OPT.is_absolute_tolerance
        tol = OPT.tol; rtol = 0;
    else
        rtol = OPT.tol; tol = 0;
    end
end
if isfield(OPT,'itmax')
    itmax = OPT.itmax;
end
if isfield(OPT,'PC')
    PC = OPT.PC;
end
if isfield(OPT,'Jwght_down') && length(OPT.Jwght_down)==nmesh
    Jwght_down = OPT.Jwght_down;
end
if isfield(OPT,'Jwght_up') && length(OPT.Jwght_up)==nmesh
    Jwght_up = OPT.Jwght_up;
end

end % END OF FUNCTION setup_mgpcg

