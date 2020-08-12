function normx = normdf(x)
% Usage: normx = normdf(x)
% 
% Purpose: Calculates norm: norm(x)/sqrt(length(x))
%
% Input:
%   x     : [vector or matrix] : vector or matrix of which norm is to be
%                                calculated
% Output:
%   normx : [scalar]           : norm of x
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH & JPM, May, 31st, 2008
%

normx = norm(double(x))/sqrt(length(x(:)));

end