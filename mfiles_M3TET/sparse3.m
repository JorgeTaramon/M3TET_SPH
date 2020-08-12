function A = sparse3(iA,jA,vA,nA,mA)
% Usage: A = sparse3(iA,jA,vA,nA,mA)
% 
% Purpose: Uses most efficient "sparse" method available (either Matlab's
%          sparse or SuiteSparse's sparse2)
%
% Input:
%   iA : [colvector] : row-index for each matrix entry
%   jA : [colvector] : col-index for each matrix entry
%   vA : [colvector] : value for each matrix entry
%   nA : [scalar]    : number of rows of assembled matrix (optional)
%   mA : [scalar]    : number of cols of assembled matrix (optional)
%
% Output:
%   A    : [sparsemat] : symmetric pos-def matrix to be factorized
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Dec 2012
% JH May 2015 :: use "exist" to see if sparse2 is available
%

if exist('sparse2','file')
    switch nargin
	case 1
	    A = sparse2(iA); % iA is a matrix when called with nargin==1
        case 3
            A = sparse2(iA,jA,vA);
        case 5
            A = sparse2(iA,jA,vA,nA,mA);
        otherwise
            error('Either 1, 3 or 5 input arguments.');
    end
else
    switch nargin
        case 1
	    A = sparse(iA); % iA is a matrix when called with nargin==1
        case 3
            A = sparse(double(iA),double(jA),vA);
        case 5
            A = sparse(double(iA),double(jA),vA,double(nA),double(mA));
        otherwise
            error('Either 1, 3 or 5 input arguments.');
    end
end

end % END OF FUNCTION sparse3
