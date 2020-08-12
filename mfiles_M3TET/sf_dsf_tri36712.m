function [N,dN] = sf_dsf_tri36712(lc,nnodel,return_format,shape_function)
% Usage: [N,dN] = sf_dsf_tri36712(lc,nnodel,return_format,shape_function)
% 
% Purpose: Calculates values of shape functions and their derivatives at
%          points with local coordinates (r,s,t=1-r-s) inside triangles.
%
% Input:
%   lc             : [matrix] : local coordinates of points (2 x npt)
%   nnodel         : [scalar] : number of nodes per triangular element
%   return_format  : [string] : 'cell' or 'matrix'
%   shape_function : [string] : 'std' or 'hier'
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
% JMT Jun 2016: added hierarchical shape functions

if nargin<3
    error('Third input argument is missing (it defines the output format)! Must be either "cell" or "matrix".');
end
if ~ismember(nnodel,[3 6 7 12])
    error(' nnodel must be 3, 6, 7, or 12');
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
                switch shape_function
                    case 'std'
                        for ip=1:npt
                            [N{ip},dN{ip}] = sf_dsf_tri6(r(ip),s(ip)); % *SUBFUNCTION*
                        end
                    case 'hier'
                        for ip=1:npt
                            [N{ip},dN{ip}] = sf_dsf_tri6h(r(ip),s(ip)); % *SUBFUNCTION*
                        end
                end

            case 7 % 7-node triangle
                for ip=1:npt
                    % VERSION (1): standard 7-node shape functions
                    [N{ip},dN{ip}] = sf_dsf_tri7(r(ip),s(ip)); % *SUBFUNCTION*
                end
                                
            case 12 % 12-node triangle
                for ip=1:npt
                    [N{ip},dN{ip}] = sf_dsf_tri12(r(ip),s(ip)); % *SUBFUNCTION*
                end
        end
        
    case 'matrix'
        dN = [];
        switch nnodel
            case 3 % 3-node triangle
                N = sf_dsf_tri3(r,s); % *SUBFUNCTION*
                
            case 6 % 6-node triangle
                switch shape_function
                    case 'std'
                        N = sf_dsf_tri6(r,s); % *SUBFUNCTION*
                    case 'hier'
                        N = sf_dsf_tri6h(r,s); % *SUBFUNCTION*
                end
                
            case 7 % 7-node triangle
                % VERSION (1): standard 7-node shape functions
                N = sf_dsf_tri7(r,s); % *SUBFUNCTION*
                
            case 12 % 12-node triangle
                N = sf_dsf_tri12(r,s); % *SUBFUNCTION*
                
        end
        
    otherwise
        error('Third input argument (defining output format) must be "cell" or "matrix".');
end

end % END OF FUNCTION sf_dsf_tri36712

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function [N,dN] = sf_dsf_tri12(r,s)
% Find shape functions and their derivatives at given points on the
% master element for a 12 node triangle

% 12-node triangle (node numbering is important)
%
%         3
%         | \
%         8  7
% s-axis  | 12 \
%         9     6
%         | 10 11 \
%         1--4--5--2
%           r-axis

t = 1-r-s;

N = [4.5 .* t .* (t-1/3) .* (t-2/3) ... % N1 at coordinate (r,s)
     4.5 .* r .* (r-1/3) .* (r-2/3) ... % N2 at coordinate (r,s)
     4.5 .* s .* (s-1/3) .* (s-2/3) ... % etc
     13.5 .* r .* t .* (t-1/3) ...
     13.5 .* r .* t .* (r-1/3) ...
	 13.5 .* s .* r .* (r-1/3) ...
	 13.5 .* s .* r .* (s-1/3) ...
	 13.5 .* t .* s .* (s-1/3) ...
	 13.5 .* t .* s .* (t-1/3) ...
     27 .* r .* s .* t.^2 ...
     27 .* r.^2 .* s .* t ...
     27 .* r .* s.^2 .* t]';

if nargout==2
    dN = [18*r + 18*s - 27*r*s - (27*r^2)/2 - (27*s^2)/2 - 11/2 ... % dN1/dr
         (27*r^2)/2 - 9*r + 1 ...                                   % dN2/dr
         0 ...                                                      % dN3/dr
         (81*r^2)/2 + 54*r*s - 45*r + (27*s^2)/2 - (45*s)/2 + 9 ...
         36*r + (9*s)/2 - 27*r*s - (81*r^2)/2 - 9/2 ...
         (9*s*(6*r - 1))/2 ...
         (27*s*(s - 1/3))/2 ...
         -(27*s*(s - 1/3))/2 ...
         (9*s*(6*r + 6*s - 5))/2 ...
         27*s*(r + s - 1)*(3*r + s - 1) ...
        -27*r*s*(3*r + 2*s - 2) ...
        -27*s^2*(2*r + s - 1); ...       
         ...
         18*r + 18*s - 27*r*s - (27*r^2)/2 - (27*s^2)/2 - 11/2 ... % dN1/ds
         0 ...                                                     % dN2/ds
         (27*s^2)/2 - 9*s + 1 ...                                  % dN3/ds
         (9*r*(6*r + 6*s - 5))/2 ...
         -(27*r*(r - 1/3))/2 ...
         (27*r*(r - 1/3))/2 ...
         (9*r*(6*s - 1))/2 ...
         (9*r)/2 + 36*s - 27*r*s - (81*s^2)/2 - 9/2 ...
         (27*r^2)/2 + 54*r*s - (45*r)/2 + (81*s^2)/2 - 45*s + 9 ...
         27*r*(r + s - 1)*(r + 3*s - 1) ...
        -27*r^2*(r + 2*s - 1) ...
        -27*r*s*(2*r + 3*s - 2)];
else
    dN = [];
end
    
end % END OF FUNCTION sf_dsf_tri12

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

function [N,dN] = sf_dsf_tri6h(r,s)
% Find hierarchical shape functions and their derivatives at given points 
% on the master element for a 6 node triangle 

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

N = [t       ... % N1 at coordinate (r,s)
     r       ... % N2 at coordinate (r,s)
     s       ... % etc
     4.*r.*t ...
     4.*r.*s ...
     4.*s.*t ]';
 
if nargout==2
    %     dN1 dN2 dN3  dN4      dN5   dN6
    dN = [-1  1   0    4*(t-r)  4*s  -4*s  ;     % w.r.t. r
          -1  0   1   -4*r      4*r   4*(t -s)]; % w.r.t. s
else
    dN = [];
end

end % END OF FUNCTION sf_dsf_tri6h

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