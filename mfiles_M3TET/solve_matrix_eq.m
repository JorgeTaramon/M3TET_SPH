function [x,profiling] = solve_matrix_eq(A,b,x0,COMM,OPTS,fidl)

% if nargin==0
%     load('StokesFlowExample','A','b','x0','Ic2f','ifree');
%     x = x0;
%     COMM.myid         = 1;
%     COMM.nsd          = 1;
%     COMM.nthreads     = 1;
%     OPTS.solver       = 'pcg';
%     OPTS.CoarseSolver = 'chol';
%     OPTS.tol          = 1e-6;
%     OPTS.is_abs_tol   = 0;
%     OPTS.itmax        = 1000;
%     OPTS.pc           = 'mg';
%     OPTS.FigNo        = 1;
% end

if size(A,1)~=size(A,2)
    error('A is not a square matrix');
else
    n = length(A);
end
if length(b)~=n
    error('Matrix A does not match rhs-vector b');
end
if nargin<6
    fidl = 1; % output on screen
end

t=tic;
switch OPTS.solver
    case 'pcg'
        if ~isfield(OPTS,'use_mutils') || OPTS.use_mutils==0
            % Have to create second half of coefficient matrix when not using MUTILS
            if max(abs(A(1,2:end)))==0
                A = A + tril(A,-1)';
            end
        end
        if ~exist('x0','var')
            x0 = zeros(size(b));
        end
        if ~strcmp(OPTS.pc,'mg')
            [x,profiling] = pcg_solver_p(A,b,x0,COMM,OPTS);
        else
            error('Multigrid preconditioner for Conjugate Gradient Solver only in viscous flow solver.');
        end
        fprintf(fidl,'(PCG,%5s,tol=%.1e)  : %7.2f sec (%4i iterations)\n',...
            profiling.pc,profiling.tol,toc(t),profiling.nit);
        
    case 'chol'
        tic
        [L,perm] = cholesky_factorization(A);
        x        = cholesky_solve(b,L,perm);
        profiling.t = toc;
        fprintf(fidl,'(chol):                %7.2f sec\n',toc(t));
        
    case 'bslash'
        tic
        if max(abs(A(1,2:end)))==0
            A = A + tril(A,-1)';
        end
        x = A \ b;
        fprintf(fidl,'(backslash):           %7.2f sec\n',toc(t));
end

end