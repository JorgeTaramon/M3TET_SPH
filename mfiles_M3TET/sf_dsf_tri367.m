function [N,dN] = sf_dsf_tri367(lc,nnodel,return_format)
% Usage: [N,dN] = sf_dsf_tri367(lc,nnodel,return_format)
% 
% Purpose: Calculates values of shape functions and their derivatives at
%          points with local coordinates (r,s,t=1-r-s) inside triangles.
%
% Input:
%   lc            : [matrix] : local coordinates of points (2 x npt)
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
% Part of M2TRI - 2D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Dec 2012
% JH Mar 2014: values returned in cells
% JH Nov 2014: function can be called in two different ways using
%    1) return_format=='cell' --> N and dN return as cell array for each 
%                                 point inside the triangle
%    2) return_format=='matrix' --> only N is returned in a matrix of size
%                                   [nnodel x npt]; dN=[]
%

if nargin<3
    error('Third input argument is missing (it defines the output format)! Must be either "cell" or "matrix".');
end
if ~ismember(nnodel,[3 6 7])
    error(' nnodel must be 3,6 or 7.');
end
r    = lc(1,:)';
s    = lc(2,:)';
npt  = length(r);
switch return_format
    case 'cell'
        N  = cell(npt,1);
        dN = cell(npt,1);
        switch nnodel
            case 3 % 3-node triangle
                for ip=1:npt
                    [N{ip},dN{ip}] = sf_dsf_tri3(r(ip),s(ip)); % *SUBFUNCTION*
                end

            case 6 % 6-node triangle
                for ip=1:npt
                    [N{ip},dN{ip}] = sf_dsf_tri6(r(ip),s(ip)); % *SUBFUNCTION*
                end

            case 7 % 7-node triangle
                for ip=1:npt
                    % VERSION (1): standard 7-node shape functions
                    [N{ip},dN{ip}] = sf_dsf_tri7(r(ip),s(ip)); % *SUBFUNCTION*

        %             % VERSION (2): standard 6-node shape functions plus central bubble fct
        %             [N{ip},dN{ip}] = sf_dsf_tri6p1(r(ip),s(ip)); % *SUBFUNCTION*
                end
        end
        
    case 'matrix'
        dN = [];
        switch nnodel
            case 3 % 3-node triangle
                N = sf_dsf_tri3(r,s); % *SUBFUNCTION*
                
            case 6 % 6-node triangle
                N = sf_dsf_tri6(r,s); % *SUBFUNCTION*
                
            case 7 % 7-node triangle
                % VERSION (1): standard 7-node shape functions
                N = sf_dsf_tri7(r,s); % *SUBFUNCTION*
                
%                 % VERSION (2): standard 6-node shape functions plus central bubble fct
%                 N = sf_dsf_tri6p1(r,s); % *SUBFUNCTION*
        end
        
    otherwise
        error('Third input argument (defining output format) must be "cell" or "matrix".');
end

end % END OF FUNCTION sf_dsf_tri367

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function [N,dN] = sf_dsf_tri7(r,s)
% Find shape functions and their derivatives at given points on the
% master element for a 7 node triangle

% 7-node triangle (node numbering is important)
%
%        3
%        | \
% s-axis 6   5
%        | 7   \
%        1 - 4 - 2
%          r-axis

t = 1-r-s;

N = [t.*(2.*t-1) + 3.*r.*s.*t ... % N1 at coordinate (r,s)
     r.*(2.*r-1) + 3.*r.*s.*t ... % N2 at coordinate (r,s)
     s.*(2.*s-1) + 3.*r.*s.*t ... % etc
     4.*r.*t - 12.*r.*s.*t ...
     4.*r.*s - 12.*r.*s.*t ...
     4.*s.*t - 12.*r.*s.*t ...
     27.*r.*s.*t           ]';

if nargout==2
    dN = [ 1-4*t+3*s*t-3*r*s ... % dN1/dr
          -1+4*r+3*s*t-3*r*s ... % dN2/dr
           3*s*t-3*r*s ...       % etc
           4*t-4*r+12*r*s-12*s*t ...
           4*s+12*r*s-12*s*t ...
          -4*s+12*r*s-12*s*t ...
          -27*r*s+27*s*t;...
          ...
           1-4*t+3*r*t-3*r*s ... % dN1/ds
           3*r*t-3*r*s ...       % dN2/ds
          -1+4*s+3*r*t-3*r*s ... % etc
          -4*r-12*r*t+12*r*s ...
           4*r-12*r*t+12*r*s ...
           4*t-4*s-12*r*t+12*r*s ...
           27*r*t-27*r*s         ];
else
    dN = [];
end
    
end % END OF FUNCTION sf_dsf_tri7

% ###################################################################

function [N,dN] = sf_dsf_tri6p1(r,s)
% Find shape functions and their derivatives at given points on the
% master element for a 7 node triangle

% 7-node triangle (node numbering is important)
%
%        3
%        | \
% s-axis 6   5
%        | 7   \
%        1 - 4 - 2
%          r-axis
% NOTE: 7th node is a bubble node. We use a STANDARD 6 node triangle 
%       augmented by a bubble centroid function whose associated variables
%       are therefore the DEVIATIONS of velocities from the quadratic 
%       solution of a 6-node triangle.
%

t = 1-r-s;

N = [t.*(2.*t-1) ... % N1 at coordinate (r,s)
     r.*(2.*r-1) ... % N2 at coordinate (r,s)
     s.*(2.*s-1) ... % etc
     4.*r.*t     ...
     4.*r.*s     ...
     4.*s.*t     ...
     27.*r.*s.*t ]';
 
if nargout==2
    %     dN1       dN2    dN3     dN4      dN5   dN6       dN7
    dN = [-(4*t-1)  4*r-1  0       4*(t-r)  4*s   -4*s      27*(s*t-r*s);  % w.r.t. r
          -(4*t-1)  0      4*s-1  -4*r      4*r    4*(t -s) 27*(r*t-r*s)]; % w.r.t. s
else
    dN = [];
end

end % END OF FUNCTION sf_dsf_tri6p1

% ###################################################################

function [N,dN] = sf_dsf_tri6(r,s)
% Find shape functions and their derivatives at given points on the
% master element for a 6 node triangle

% 6-node triangle (node numbering is important)
%
%        3
%        | \
% s-axis 6   5
%        |     \
%        1 - 4 - 2
%          r-axis
%
t = 1-r-s;

N = [t.*(2.*t-1) ... % N1 at coordinate (r,s)
     r.*(2.*r-1) ... % N2 at coordinate (r,s)
     s.*(2.*s-1) ... % etc
     4.*r.*t     ...
     4.*r.*s     ...
     4.*s.*t    ]';
 
if nargout==2
    %     dN1       dN2    dN3    dN4       dN5    dN6
    dN = [-(4*t-1)  4*r-1  0       4*(t-r)  4*s   -4*s  ;    % w.r.t. r
          -(4*t-1)  0      4*s-1  -4*r      4*r   4*(t -s)]; % w.r.t. s
else
    dN = [];
end

end % END OF FUNCTION sf_dsf_tri6

% ###################################################################

function [N,dN] = sf_dsf_tri3(r,s)
% Find shape functions and their derivatives at given points on the
% master element for 3 node triangle

% 3-node triangle (node numbering is important)
%
%        3
%        | \
% s-axis |   \
%        |     \
%        1 - - - 2     
%          r axis -->
%

t  = 1-r-s;

N  = [t ... % N1 at coordinate (r,s)
      r ... % N2 at coordinate (r,s)
      s]';   % N3 at coordinate (r,s)

if nargout==2
    %     dN1 dN2 dN3
    dN = [-1   1   0;  % w.r.t. r
          -1   0   1]; % w.r.t. s
else
    dN = [];
end

end % END OF FUNCTION sf_dsf_tri3