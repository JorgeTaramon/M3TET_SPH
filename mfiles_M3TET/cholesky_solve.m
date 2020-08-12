function x = cholesky_solve(rhs,L,perm)
% Usage: x = cholesky_solve(rhs,L,perm)
% 
% Purpose: Uses most efficient method available to perform forward-backward
%          substitution (part 2 of Cholesky direct solver)
%
% Input:
%   rhs  : [colvector] : right-hand-side of A*x = rhs
%   L    : [sparsemat] : factorized matrix (L*L'=A)
%   perm : [rowvector] : vector of row/column permutations during
%                        factorization
%
% Output:
%   x    : [colvector] : solution vector of A*x = rhs
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Dec 2012
% JH May 2015 :: use "exist" to see if SuiteSparse is available
%

if exist('cs_transpose','file') && ...
   exist('cs_symperm','file')
    % Use SuiteSparse
    if nargin==3
        x       = zeros(size(rhs));
        x(perm) = cs_ltsolve(L,cs_lsolve(L,rhs(perm)));
    else
        x = cs_ltsolve(L,cs_lsolve(L,rhs));
    end
else
    if nargin==3
        x       = zeros(size(rhs));
        x(perm) = L'\ (L \rhs(perm));
    else
        x = L'\ (L \rhs);
    end
end

end % END OF FUNCTION cholesky_solve