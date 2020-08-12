function [L,perm] = cholesky_factorization(A)
% Usage: [L,perm] = cholesky_factorization(A)
% 
% Purpose: Uses most efficient method available to factorize sparse matrix
%          (part 1 of Cholesky direct solver)
%
% Input:
%   A    : [sparsemat] : symmetric pos-def matrix to be factorized
%
% Output:
%   L    : [sparsemat] : factorized matrix (L*L'=A)
%   perm : [rowvector] : vector of row/column permutations during
%                        factorization
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

if exist('lchol','file')
    % Use SuiteSparse
    perm = amd(A);
    A    = cs_transpose(A);
    A    = cs_symperm(A,perm);
    A    = cs_transpose(A);
    L    = lchol(A);
else
    % Use Matlab's chol
    [L,~,perm] = chol(tril(A),'lower','vector');    
end

end % END OF FUNCTION cholesky_factorization