function [x,profiling] = pcg_solver_p(A,b,x,COMM,OPTS)
%
% Purpose: Conjugate Gradient solution algorithm to solve the 
%          symmetric and pos-def matrix equation A*x = b; parallel version
% Input:
%   A    : [sparsemat] : symmetric pos-def coefficient matrix
%   b    : [colvector] : right-hand-side vector of A*x = b
%   x    : [colvector] : guess for solution vector
%   COMM : [structure] : communication data
%   OPTS : [structure] : options for conjugate gradient solver
%
% Output:
%   x         : [colvector] : solution x
%   profiling : [structure] : profiling information
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)

% JH Apr 2011
% JH Dec 2014 : added MUTILS' spmv
%

check_input_args(); % *NESTED FUNCTION

if nargin<4; OPTS=[]; end
OPTS = setup_pcg(OPTS); % *SUBFUNCTION*

nsd  = COMM.nsd;
if nsd>1
    NB           = COMM.NB;     % list of neighboring SDs
    sum_SDB      = COMM.sum_SDB; % sum values at SD boundaries (function handle)
    mynod_SDB    = COMM.mynod_SDB_f;   % subdomain boundary information for all neighbors
    unique_nodes = COMM.unique_nodes_f; % unique list of domain nodes
end

t=tic;

% =========================================================================
% PREPARE PRECONITIONER
% 'diag' : diagonal scaling (Jacobi preconditioner)
% 'noPC' : no preconditioning (CG algorithm)
% 'symGS': symmetric Gauss-Seidel
% '0chol': zero-infill incomplete Cholesky factorization
% 'ichol': drop level incomplete Cholesky factorization
switch OPTS.pc
    case 'diag'
        AD  = full(diag(A));
        if nsd>1; AD = sum_SDB(AD,NB,mynod_SDB); end
        ADi = 1./AD;
    case 'noPC'
        % nothing to prepare
    case 'symGS'
        if nsd>1; error('pc==symGS not supported in parallel mode.'); end
        ADL = tril(A);
        ADU = ADL';
        AD  = full(diag(A));
    otherwise
        if nsd>1; error('pc==0chol & ichol not supported in parallel mode.'); end
        switch OPTS.pc
            case '0chol'
                OPT_ichol.type    = 'nofill';
            case 'ichol'
                OPT_ichol.type    = 'ict';
                OPT_ichol.droptol = OPTS.dropfactor_ichol * max(max(abs(A)));
        end
        perm = amd(A);
        L    = cs_transpose(A);
        L    = cs_symperm(L,perm);
        L    = cs_transpose(L);
        L    = ichol(L,OPT_ichol);
end
% =========================================================================

if OPTS.use_mutils && exist('sparse_convert','file')
    % If MUTILS is available: convert the sparse matrix to the more efficient
    % format (see help sparse_convert)
    opts_spc.symmetric = 1;
    opts_spc.nthreads  = COMM.nthreads;
    A                  = sparse_convert(tril(A),opts_spc);
end

% Calculate initial residual
% ==========================
r      = b - mat_x_vec(A,x); % *SUBFUNCTION*
if nsd>1; r = sum_SDB(r,NB,mynod_SDB); end
norm_r = normdf_p(r);  % norm of residual vector

if OPTS.is_abs_tol
    tol = OPTS.tol;       % tolerance for convergence
else
    tol = OPTS.tol*norm_r;  % tolerance for convergence
end

if norm_r<tol
    profiling.nit    = 0;
    profiling.norm_r = norm_r;
    profiling.time   = toc(t);
    profiling.pc     = OPTS.pc;
    profiling.tol    = tol;
    return
end
    
norm_r_vec    = zeros(OPTS.itmax+1,1);
norm_r_vec(1) = norm_r;
if OPTS.norm_res_vs_time
    time    = zeros(OPTS.itmax+1,1);
    time(1) = toc(t);
end
if OPTS.FigNo && nsd==1
    figure(OPTS.FigNo);ph=[];clf;
    if OPTS.norm_res_vs_time
        semilogy([0 time(1)],[norm_r_vec(1) norm_r_vec(1)],'ro-'); set(gca,'FontSize',14);
        xlabel('Time (sec)'); ylabel('Norm of residual vector');
    else
        ph=semilogy(norm_r_vec(1),'k.-'); set(gca,'FontSize',14);
        xlabel('Iteration'); ylabel('Norm of residual vector');
    end
    set(gca,'YLim',[10^floor(log10(tol)) 10^ceil(log10(norm_r))]); hold all
    semilogy([0 50],[tol tol],'r--');
end

% Begin of CG iterations
% ======================
for it=1:OPTS.itmax
    
    % Preconditioning (approximate error of current x)
    % ================================================
    d = precondition; % *NESTED FUNCTION*
    
    % get new search direction q
    % ==========================
    if it==1
        % Initialize Conjugate Gradients
        % ==============================
        q  = d; % define FIRST search direction q

    else
        rdlast = rd; % save old rd (for calculating beta)
        % Make new search direction q A-orthogonal to all previous q's
        beta = dot_p(r-rlast,d)/rdlast;  % Polak-Ribiere version 1
%         beta = dot_p(d-dlast,r)/rdlast;  % Polak-Ribiere version 2
        q    = d + beta*q; % calculate NEW search direction
    end

    rd    = dot_p(r,d);  % numerator for calculating alpha
    Aq    = mat_x_vec(A,q); % *SUBFUNCTION*
    if nsd>1; Aq = sum_SDB(Aq,NB,mynod_SDB); end
    qAq   = dot_p(q,Aq); % denominator in calculating alpha
    
    rlast = r;         % needed for Polak-Ribiere version 1
%     dlast = d;         % needed for Polak-Ribiere version 2
    alpha = rd/qAq;    % step size in direction q

    % Update solution and residual
    % ============================
    x     = x + alpha * q;  % update solution
    r     = r - alpha * Aq; % update residual
    
    % Residual norm and plotting
    % ==========================
    norm_r           = normdf_p(r); % norm of residual vector
    norm_r_vec(it+1) = norm_r;
    if OPTS.norm_res_vs_time; time(it+1) = toc(t); end
    if OPTS.FigNo && nsd==1
        figure(OPTS.FigNo);if ~isempty(ph);delete(ph);end
        if OPTS.norm_res_vs_time
            ph=semilogy(time(1:it+1),norm_r_vec(1:it+1),'k.-');
        else
            ph=semilogy(norm_r_vec(1:it+1),'k.-');
        end
        drawnow
    end
    
    % Check if solution converged
    % ===========================
    if norm_r<tol
        break
    end
end % end of CG loop

profiling.nit  = it;
profiling.norm_r = norm_r_vec(1:it+1);
if OPTS.norm_res_vs_time
    profiling.time = time(1:it+1);
else
    profiling.time = toc(t);
end
profiling.pc  = OPTS.pc;
profiling.tol = tol;

if OPTS.FigNo
    figure(OPTS.FigNo);
    if nsd>1
        set(gcf,'Visible','off');
        if OPTS.norm_res_vs_time
            semilogy(time(1:it+1),norm_r_vec(1:it+1),'k.-'); set(gca,'FontSize',14);
            xlabel('Iteration'); ylabel('Norm of residual vector');
        else
            semilogy(norm_r_vec(1:it+1),'k.-'); set(gca,'FontSize',14);
            xlabel('Iteration'); ylabel('Norm of residual vector');
        end
        set(gca,'YLim',[tol 10^ceil(log10(max(norm_r_vec)))]); hold all
        drawnow
    end
    if OPTS.norm_res_vs_time
        text(time(it+1),8*norm_r_vec(it+1),['pc=' OPTS.pc]);
        text(time(it+1),2*norm_r_vec(it+1),['nit=' num2str(it+1)]);
    else
        text(it+1,8*norm_r_vec(it+1),['pc=' OPTS.pc]);
        text(it+1,2*norm_r_vec(it+1),sprintf('t=%0.2e sec',profiling.time));
    end
    if nsd>1
        print('-dpng','-r150',[COMM.prefix '_PCGconvergence']);
    end
end

% =========================================================================
%                           NESTED FUNCTIONS
% =========================================================================

function rms = normdf_p(a)
    if nsd>1
        rms = COMM.normdf(a,unique_nodes);
    else
        rms = normdf(a);
    end
end % END OF NESTED FUNCTION normdf_p

% =========================================================================

function ab = dot_p(a,b)
    if nsd>1
        ab = COMM.dot(a,b,unique_nodes);
    else
        ab = dot(a,b);
    end
end % END OF NESTED FUNCTION dot_p

% =========================================================================

function d = precondition
    switch OPTS.pc
        case 'diag' % Diagonal scaling (Jacobi)
            d = r.*ADi;
        
        case 'noPC' % No preconditioning (i.e. CG algorithm)
            d = r;

        case 'symGS' % Symmetric Gauss-Seidel
            w = ADL \ r; % FORWARD SOLVE
            w = w .* AD; % multiply diagonal
            d = ADU \ w; % BACKWARD SOLVE
        
        otherwise
            d       = zeros(size(r));
            d(perm) = cs_ltsolve(L,cs_lsolve(L,r(perm)));
%             d = L' \ (L\ r); % SLOWER        
    end
end % END OF NESTED FUNCTION precondition

% =========================================================================

function check_input_args()
    % check input
    n = length(b);
    if size(x)~=size(b)
        error('Initial guess vector does not match rhs vector.');
    end
    if n~=size(A,1) || n~=size(A,2)
        error('Length of RHS vector does not match coefficient matrix');
    end
end % END OF NESTED FUNCTION check_input_args

end % END OF FUNCTION pcg_solver_p

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function OPTS = setup_pcg(OPTS)

if ~isfield(OPTS,'tol')
    OPTS.tol = 1e-8; % relative or absolute tolerance (defined below)
end

if ~isfield(OPTS,'is_abs_tol')
    OPTS.is_abs_tol = 1; % absolute or relative tolerance
end

if ~isfield(OPTS,'itmax')
    OPTS.itmax = 300;   % maximum number of iterations
end

if ~isfield(OPTS,'FigNo')
    OPTS.FigNo = 0;     % set =0 to disable plotting
end

if ~isfield(OPTS,'norm_res_vs_time')
    OPTS.norm_res_vs_time = 0; % =0 --> plot norm residual vs # iteration
                               % =1 --> plot norm residual vs computer time
end

if isfield(OPTS,'pc')
    if ~ismember(OPTS.pc,{'0chol','ichol','symGS','diag','noPC'})
        error(' %s is not a supported preconditioner in "pcg_solver_p"');
    end
else
    % '0chol': zero-infill incomplete Cholesky factorization
    % 'ichol': drop level incomplete Cholesky factorization
    % 'symGS': symmetric Gauss-Seidel
    % 'diag' : diagonal scaling (Jacobi preconditioner)
    % 'noPC' : no preconditioning (CG algorithm)
    OPTS.pc = 'diag';
end

if strcmp(OPTS.pc,'ichol') && ~isfield(OPTS,'dropfactor_ichol')
    OPTS.dropfactor_ichol = 5e-8;
end

if ~isfield(OPTS,'use_mutils')
    OPTS.use_mutils = 0;
end

end % END OF SUBFUNCTION setup_pcg

% #########################################################################

function x = mat_x_vec(A,b)

% Times in brackets are for 65 muliplications with nU==225370
if isstruct(A)
    x = spmv(A,b); % performs: 4.1-4.5 sec
else
    if ~any(A(1,2:length(A)))
        x = A*b + A'*b - full(diag(A)).*b; % performs: 10.1-10.6 sec
%         x = A*b + tril(A,-1)'*b; % performs: 17.5-17.8 sec
    else
        x = A*b; % performs: 6.4-6.9 sec
    end
end

end % END OF SUBFUNCTION mat_x_vec