function [N,dN] = sf_dsf_tet410(lc,nnodel,return_format)
% Usage: [N,dN] = sf_dsf_tet410(lc,nnodel,return_format)
% 
% Purpose: Calculates values of shape functions and their derivatives at
%          points with local coordinates (r,s,t,u=1-r-s-t) inside tetrahedra.
%
% Input:
%   lc            : [matrix] : local coordinates of points (3 x npt)
%   nnodel        : [scalar] : number of nodes per triangular element
%   return_format : [string] : 'cell' or 'matrix'
%
% Output (if "return_format"=='cell')
%   N      : [cell]   : shape fct values (npt x 1)
%   dN     : [cell]   : shape fct derivatives (npt x 1)
%
% Output (if "return_format"=='matrix')
%   N      : [matrix] : shape fct values (npt x nnodel)
%   dN     : []       : NOT CALCULATED, RETURNED EMPTY
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Dec 2012
% JH Apr 2014: values returned in cells
% JH Nov 2014: function can be called in two different ways using
%    1) return_format=='cell' --> N and dN return as cell array for each 
%                                 point inside the triangle
%    2) return_format=='matrix' --> only N is returned in a matrix of size
%                                   [nnodel x npt]; dN=[]
%

error('use function sf_dsf_tet.m')

if nargin<3
    error('Third input argument is missing (it defines the output format)! Must be either "cell" or "matrix".');
end
if ~ismember(nnodel,[4 10])
    error(' nnodel must be 4 or 10.');
end
r    = lc(1,:)';
s    = lc(2,:)';
t    = lc(3,:)';
npt  = length(r);
switch return_format
    case 'cell'
        N  = cell(npt,1);
        dN = cell(npt,1);
        switch nnodel
            case 4 % 4-node tetrahedron
                for ip=1:npt
                    [N{ip},dN{ip}] = sf_dsf_tet4(r(ip),s(ip),t(ip)); % *SUBFUNCTION*
                end

            case 10 % 10-node tetrahedron
                for ip=1:npt
                    [N{ip},dN{ip}] = sf_dsf_tet10(r(ip),s(ip),t(ip)); % *SUBFUNCTION*
                end
        end
    case 'matrix'
        dN = [];
        switch nnodel
            case 4 % 4-node tetrahedron
                N = sf_dsf_tet4(r,s,t); % *SUBFUNCTION*
            case 10 % 10-node tetrahedron
                N = sf_dsf_tet10(r,s,t); % *SUBFUNCTION*
        end
end

end % END OF FUNCTION sf_dsf_tet410

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function [N,dN] = sf_dsf_tet10(r,s,t)
% Find shape functions and their derivatives at given points on the
% master element for a 10 node tetrahedron

% Node notation taken from Hughes' book, p.171;
% see also Zienkiewicz (4th) Volume 1, p.137
% local coords of the 10 nodes:
% Node  r  s  t
%  1    1  0  0
%  2    0  1  0
%  3    0  0  1
%  4    0  0  0
%  5   .5 .5  0
%  6    0 .5 .5
%  7    0  0 .5
%  8   .5  0  0
%  9   .5  0 .5
% 10    0 .5  0

r   = r(:);
s   = s(:);
t   = t(:);
u   = 1-r-s-t;


% SHAPE FUNCTION VALUES
N  = [r.*(2.*r-1) ... % 4 vertex nodes
      s.*(2.*s-1) ...
      t.*(2.*t-1) ...
      u.*(2.*u-1) ...
      4.*r.*s ... % 6 edge nodes
      4.*s.*t ...
      4.*t.*u ...
      4.*r.*u ...
      4.*r.*t ...
      4.*s.*u]';

% DERIVATIVES
if nargout==2
    dN      = zeros(10,3); % derivatives (3 for each node)

    dN(1,1) = 4*r-1; % dN1/dr
    dN(1,2) = 0;     % dN1/ds
    dN(1,3) = 0;     % dN1/dt

    dN(2,1) = 0;     % dN2/dr
    dN(2,2) = 4*s-1; % dN2/ds
    dN(2,3) = 0;     % dN2/dt

    dN(3,1) = 0;     % dN3/dr
    dN(3,2) = 0;     % dN3/ds
    dN(3,3) = 4*t-1; % dN3/dt

    dN(4,1) = -3+4*(1-u); % dN/dr
    dN(4,2) = -3+4*(1-u); % dN/ds
    dN(4,3) = -3+4*(1-u); % dN/dt

    dN(5,1) = 4*s;   % dN5/dr
    dN(5,2) = 4*r;   % dN5/ds
    dN(5,3) = 0;     % dN5/dt

    dN(6,1) = 0;     % dN6/dr
    dN(6,2) = 4*t;   % dN6/ds
    dN(6,3) = 4*s;   % dN6/dt

    dN(7,1) = -4*t;    % dN7/dr
    dN(7,2) = -4*t;    % dN7/ds
    dN(7,3) = 4*(u-t); % dN7/dt

    dN(8,1) = 4*(u-r); % dN8/dr
    dN(8,2) = -4*r;    % dN8/ds
    dN(8,3) = -4*r;    % dN8/dt

    dN(9,1) = 4*t;   % dN9/dr
    dN(9,2) = 0;     % dN9/ds
    dN(9,3) = 4*r;   % dN9/dt

    dN(10,1) = -4*s;    % dN/dr
    dN(10,2) = 4*(u-s); % dN/ds
    dN(10,3) = -4*s;    % dN/dt
end

end % END OF FUNCTION sf_dsf_tet10

% #########################################################################

function [N,dN] = sf_dsf_tet4(r,s,t)
% Find shape functions and their derivatives at given points on the
% master element for 4 node tetrahedron

% Node notation taken from Hughes' book, p.170
% local coords of the 4 nodes:
% Node  r  s  t
%  1    1  0  0
%  2    0  1  0
%  3    0  0  1
%  4    0  0  0

r   = r(:);
s   = s(:);
t   = t(:);

% SHAPE FUNCTION VALUES
N  = [r s t 1-r-s-t]';

% DERIVATIVES
if nargout==2
    dN      = zeros(4,3); % derivatives (3 for each node)

    dN(1,1) = 1; % dN1/dr
    dN(1,2) = 0; % dN1/ds
    dN(1,3) = 0; % dN1/dt

    dN(2,1) = 0; % dN2/dr
    dN(2,2) = 1; % dN2/ds
    dN(2,3) = 0; % dN2/dt

    dN(3,1) = 0; % dN3/dr
    dN(3,2) = 0; % dN3/ds
    dN(3,3) = 1; % dN3/dt

    dN(4,1) = -1; % dN/dr
    dN(4,2) = -1; % dN/ds
    dN(4,3) = -1; % dN/dt
end

end % END OF FUNCTION sf_dsf_tet4