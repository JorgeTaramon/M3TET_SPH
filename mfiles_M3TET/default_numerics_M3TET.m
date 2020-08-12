function SETTINGS = default_numerics_M3TET()
% Usage: SETTINGS = default_numerics_M3TET()
% 
% Purpose: Defines some default settings for mesh generation and thermal
%          advection/diffusion.
%
% Input:
%   none
%
% Output:
%   SETTINGS : [structure] : model parameters
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Jan 2012
% JH Feb 2016 : updated fieldnames
%


%==========================================================================
% DEFAULT CODE SETUP
%==========================================================================

SETTINGS.disp_profiling = 1;

%==========================================================================
% DOMAIN BOUNDARY INDICES
SETTINGS.DB_indices    = cell(2,1);
SETTINGS.DB_indices{1} = 301;
SETTINGS.DB_indices{2} = 306;
  % Name the boundary indices of corners (101:199), edges (201:299) and
  % surfaces (301:399) for each domain faces
  % A sphere, however, only has a surface at bottom (301) and top (306)
%==========================================================================


%==========================================================================
% OPTIONS FOR "thermal3d_p"

% SETTINGS FOR ELEMENT ASSEMBLY
OPTS_T.eltype = 'quadratic'; % linear or quadratic
OPTS_T.lumping  = 'no';    % 'yes' or 'no'
OPTS_T.nip    =     4; % Number of integration points (1, 4, 5, 8 or 11)
OPTS_T.nelblk =  1600; % Number elements per integration-loop

% SOLVER FOR MATRIX EQUATIONS
% Specify numerical solution algorithm for the 3 matrix eq to be solved
%   'chol'   :: Cholesky direct solver (SuiteSparse if available)
%   'bslash' :: Matlab's backslash
%   'pcg'    :: Preconditioned Conjugate Gradient solver
OPTS_T.solver = 'pcg'; % 'pcg', 'chol', or 'bslash'
if strcmp(OPTS_T.solver,'pcg')
    % IF solver=='pcg' you need to define a preconditioning method:
    %   pc=='none'  :: no preconditioning (d=r)
    %   pc=='diag'  :: Jacobi diagonal scaling
    %   pc=='symGS' :: symmetric Gauss-Seidel
    %   pc=='ichol' :: droplevel incomplete Cholesky factorization
    %   pc=='0chol' :: zero-infill incomplete Cholesky factorization
    %   pc=='mg'    :: single V-cycle of multigrid
    OPTS_T.pc    = 'diag'; % Define preconditioner for CG
    if strcmp(OPTS_T.pc,'ichol')
        OPTS_T.dropfactor_ichol = 5e-8;
    end
    OPTS_T.tol   = 1e-8; % tolerance for solution
    OPTS_T.is_abs_tol = 0; % 0--> tol relative to initial residual (rhs - K*T0)
                           % 1--> tol is absolute (dangerous!!!!)
    OPTS_T.itmax = 500; % maximum number of CG iterations
    OPTS_T.FigNo = 0;   % figure number for convergence plot (0 disables plot)
    OPTS_T.norm_res_vs_time = 0;
        % 0--> plot is norm(residual) vs. iteration
        % 1--> plot is norm(residual) vs. time
    OPTS_T.use_mutils       = 1; 
        % 1--> use mutils for sparse-matrix-vector-multiplications
        % 0--> use Matlab (A*b)
end

% DEFINE HOW THE ELEMENT PROPERTIES ARE EVALUATED
OPTS_T.method_eval_dens = 'interp_nodal';
OPTS_T.method_eval_cond = 'elem_phases';
OPTS_T.method_eval_Cp   = 'elem_phases';
OPTS_T.method_eval_dQdt = 'zero_dQdt';
SETTINGS.OPTS_T = OPTS_T;
%==========================================================================


%==========================================================================
% OPTIONS FOR "diffusion3d_p"

% SETTINGS FOR ELEMENT ASSEMBLY
OPTS_D.eltype   = 'linear'; % linear or quadratic
OPTS_D.lumping  = 'yes';    % 'yes' or 'no'
OPTS_D.nip      = 4; % Number of integration points (1, 4, 5, 8 or 11)
OPTS_D.nelblk   = 1600; % Number elements per integration-loop
% SOLVER FOR MATRIX EQUATIONS
% Specify numerical solution algorithm for the 3 matrix eq to be solved
%   'chol'   :: Cholesky direct solver (SuiteSparse if available)
%   'bslash' :: Matlab's backslash
%   'pcg'    :: Preconditioned Conjugate Gradient solver
OPTS_D.solver = 'pcg'; % 'pcg', 'chol', or 'bslash'
if strcmp(OPTS_D.solver,'pcg')
    % IF solver=='pcg' you need to define a preconditioning method:
    %   pc=='none'  :: no preconditioning (d=r)
    %   pc=='diag'  :: Jacobi diagonal scaling
    %   pc=='symGS' :: symmetric Gauss-Seidel
    %   pc=='ichol' :: droplevel incomplete Cholesky factorization
    %   pc=='0chol' :: zero-infill incomplete Cholesky factorization
    %   pc=='mg'    :: single V-cycle of multigrid
    OPTS_D.pc    = 'diag'; % Define preconditioner for CG
    if strcmp(OPTS_D.pc,'ichol')
        OPTS_D.dropfactor_ichol = 5e-8;
    end
    OPTS_D.tol   = 1e-8; % tolerance for solution
    OPTS_D.is_abs_tol = 0; % 0--> tol relative to initial residual (rhs - K*T0)
                           % 1--> tol is absolute (dangerous!!!!)
    OPTS_D.itmax = 500; % maximum number of CG iterations
    OPTS_D.FigNo = 0;   % figure number for convergence plot (0 disables plot)
    OPTS_D.norm_res_vs_time = 0;
        % 0--> plot is norm(residual) vs. iteration
        % 1--> plot is norm(residual) vs. time
    OPTS_D.use_mutils       = 1; 
        % 1--> use mutils for sparse-matrix-vector-multiplications
        % 0--> use Matlab (A*b)
end
% Time approximation scheme (theta-time stepping rule)
% theta == 0   :: forward difference (EXPLICIT)
%       == 1   :: backward difference (FULL IMPLICIT)
%       == 1/2 :: midpoint rule; Crank-Nicholson (IMPLICIT)
%       == 2/3 :: Galerkin (IMPLICIT)
OPTS_D.theta    = 1; % time approximation scheme
SETTINGS.OPTS_D = OPTS_D;
%==========================================================================


%==========================================================================
% THERMAL (AND COMPOSITIONAL) ADVECTION
OPTS_SLM.method_BT = 'PC2';
    % method for evaluating origin point during back tracking
    % 'EU1' == Euler (very inaccurate !!!!)
    % 'PC2' == Predictor-Corrector
    % 'RK4' == Runge-Kutta 4th order
OPTS_SLM.method_interp = 'quadratic';
    % method for interpolating temperature at origin point
    % 'linear', 'quadratic' or 'cubic'
OPTS_SLM.monotonic   = 1;
    % only if method_interp=='quadratic' or 'cubic':
    % 1 == cutoff interpolated values so that they do not exceed
    %      the surrounding nodal values
OPTS_SLM.method_wght = 'dist';
    % Weighting for cubic interpolation
    % 'dist' or 'ones'
OPTS_SLM.nelblk   = 10000;
    % Number of elements processed at once in cubic interpolation
SETTINGS.OPTS_SLM = OPTS_SLM;
%==========================================================================

end