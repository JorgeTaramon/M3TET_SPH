function [N,dN] = sf_dsf_tet(lc,nnodel,return_format)
% Usage: [N,dN] = sf_dsf_tet(lc,nnodel,return_format)
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
% JH Mar 2016: added 20- and 21-node tetrahedra
% JH Jun 2017: added 14- and modified 15-node tetrahedra
% JMT Jul 2017: fixed bug in SFs for 14-node tetrahedra. I suspect that 
%               SFs for 15-node tetrahedra also need to be debugged
%

if nargin<3
    error('Third input argument is missing (it defines the output format)! Must be either "cell" or "matrix".');
end
if ~ismember(nnodel,[4 10 11 14 15 20 21])
    error(' nnodel must be 4, 10, 11, 14, 15, 20 or 21.');
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
            case 11 % 11-node tetrahedron
                for ip=1:npt
                    [N{ip},dN{ip}] = sf_dsf_tet11(r(ip),s(ip),t(ip)); % *SUBFUNCTION*
                end
            case 14 % 14-node tetrahedron
                for ip=1:npt
                    [N{ip},dN{ip}] = sf_dsf_tet14(r(ip),s(ip),t(ip)); % *SUBFUNCTION*
                end
            case 15 % 15-node tetrahedron
                warning('This function might need to be debugged');
                for ip=1:npt
                    [N{ip},dN{ip}] = sf_dsf_tet15(r(ip),s(ip),t(ip)); % *SUBFUNCTION*
                end
            case 20 % 20-node tetrahedron
                for ip=1:npt
                    [N{ip},dN{ip}] = sf_dsf_tet20(r(ip),s(ip),t(ip)); % *SUBFUNCTION*
                end
            case 21 % 21-node tetrahedron
                for ip=1:npt
                    [N{ip},dN{ip}] = sf_dsf_tet21(r(ip),s(ip),t(ip)); % *SUBFUNCTION*
                end
        end
    case 'matrix'
        dN = [];
        switch nnodel
            case 4 % 4-node tetrahedron
                N = sf_dsf_tet4(r,s,t); % *SUBFUNCTION*
            case 10 % 10-node tetrahedron
                N = sf_dsf_tet10(r,s,t); % *SUBFUNCTION*
            case 11 % 11-node tetrahedron
                N = sf_dsf_tet11(r,s,t); % *SUBFUNCTION*
            case 14 % 14-node tetrahedron
                N = sf_dsf_tet14(r,s,t); % *SUBFUNCTION*
            case 15 % 15-node tetrahedron
                N = sf_dsf_tet15(r,s,t); % *SUBFUNCTION*
            case 20 % 20-node tetrahedron
                N = sf_dsf_tet20(r,s,t); % *SUBFUNCTION*
            case 21 % 21-node tetrahedron
                N = sf_dsf_tet21(r,s,t); % *SUBFUNCTION*
        end
end

end % END OF FUNCTION sf_dsf_tet

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function [N,dN] = sf_dsf_tet21(r,s,t)
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
%  top be checked
% 
r   = r(:)';
s   = s(:)';
t   = t(:)';
u   = 1-r-s-t;

% SHAPE FUNCTION VALUES  
N       = zeros(21,length(r));
N( 1,:) = r.*(3.*r - 2).*((3.*r)./2 - 1/2) - (64.*r.*s.*t.*u)./5; % 4 vertex nodes
N( 2,:) = s.*(3.*s - 2).*((3.*s)./2 - 1/2) - (64.*r.*s.*t.*u)./5;
N( 3,:) = t.*(3.*t - 2).*((3.*t)./2 - 1/2) - (64.*r.*s.*t.*u)./5;
N( 4,:) = u.*(3.*u - 2).*((3.*u)./2 - 1/2) - (64.*r.*s.*t.*u)./5;
N( 5,:) = r.*s.*((27.*r)./2 - 9/2) - (64.*r.*s.*t.*u)./5; % 12 edge nodes
N( 6,:) = r.*s.*((27.*s)./2 - 9/2) - (64.*r.*s.*t.*u)./5;
N( 7,:) = s.*t.*((27.*s)./2 - 9/2) - (64.*r.*s.*t.*u)./5;
N( 8,:) = s.*t.*((27.*t)./2 - 9/2) - (64.*r.*s.*t.*u)./5;
N( 9,:) = t.*u.*((27.*t)./2 - 9/2) - (64.*r.*s.*t.*u)./5;
N(10,:) = t.*u.*((27.*u)./2 - 9/2) - (64.*r.*s.*t.*u)./5;
N(11,:) = r.*u.*((27.*u)./2 - 9/2) - (64.*r.*s.*t.*u)./5;
N(12,:) = r.*u.*((27.*r)./2 - 9/2) - (64.*r.*s.*t.*u)./5;
N(13,:) = r.*t.*((27.*r)./2 - 9/2) - (64.*r.*s.*t.*u)./5;
N(14,:) = r.*t.*((27.*t)./2 - 9/2) - (64.*r.*s.*t.*u)./5;
N(15,:) = s.*u.*((27.*u)./2 - 9/2) - (64.*r.*s.*t.*u)./5;
N(16,:) = s.*u.*((27.*s)./2 - 9/2) - (64.*r.*s.*t.*u)./5;
N(17,:) = 27.*s.*t.*u - (64.*r.*s.*t.*u)./5; % 4 face nodes
N(18,:) = 27.*r.*t.*u - (64.*r.*s.*t.*u)./5;
N(19,:) = 27.*r.*s.*u - (64.*r.*s.*t.*u)./5;
N(20,:) = 27.*r.*s.*t - (64.*r.*s.*t.*u)./5;
N(21,:) = 256.*r.*s.*t.*u; % central bubble node

% DERIVATIVES
if nargout==2
    dN       = zeros(21,3); % derivatives (3 for each node)
    
    dN( 1,:) = [(27*r^2)/2 + (128*r*s*t)/5 - 9*r + (64*s^2*t)/5 + (64*s*t^2)/5 - (64*s*t)/5 + 1 ...
                (64*r*t*(s - u))/5 ...
                (64*r*s*(t - u))/5];

    dN( 2,:) = [(64*s*t*(r - u))/5 ...
                (64*r^2*t)/5 + (128*r*s*t)/5 + (64*r*t^2)/5 - (64*r*t)/5 + (27*s^2)/2 - 9*s + 1 ...
                (64*r*s*(t - u))/5];

    dN( 3,:) = [(64*s*t*(r - u))/5 ...
                (64*r*t*(s - u))/5 ...
                (64*r^2*s)/5 + (64*r*s^2)/5 + (128*r*s*t)/5 - (64*r*s)/5 + (27*t^2)/2 - 9*t + 1];

    dN( 4,:) = [(64*s*t^2)/5 - 27*r*s - 27*r*t - (199*s*t)/5 - 18*u + (64*s^2*t)/5 - (27*r^2)/2 - (27*s^2)/2 - (27*t^2)/2 + (128*r*s*t)/5 + 25/2 ...
                (64*r*t^2)/5 - 27*r*s - (199*r*t)/5 - 27*s*t - 18*u + (64*r^2*t)/5 - (27*r^2)/2 - (27*s^2)/2 - (27*t^2)/2 + (128*r*s*t)/5 + 25/2 ...
                (64*r*s^2)/5 - (199*r*s)/5 - 27*r*t - 27*s*t - 18*u + (64*r^2*s)/5 - (27*r^2)/2 - (27*s^2)/2 - (27*t^2)/2 + (128*r*s*t)/5 + 25/2];

    dN( 5,:) = [(s*(270*r - 128*t + 256*r*t + 128*s*t + 128*t^2 - 45))/10 ...
                (r*(135*r + 128*s*t - 128*t*u - 45))/10 ...
                (64*r*s*(t - u))/5];

    dN( 6,:) = [(s*(135*s + 128*r*t - 128*t*u - 45))/10 ...
                (r*(270*s - 128*t + 128*r*t + 256*s*t + 128*t^2 - 45))/10 ...
                (64*r*s*(t - u))/5];

    dN( 7,:) = [(64*s*t*(r - u))/5 ...
                (t*(270*s - 128*r + 256*r*s + 128*r*t + 128*r^2 - 45))/10 ...
                (s*(135*s + 128*r*t - 128*r*u - 45))/10];

    dN( 8,:) = [(64*s*t*(r - u))/5 ...
                (t*(135*t + 128*r*s - 128*r*u - 45))/10 ...
                (s*(270*t - 128*r + 128*r*s + 256*r*t + 128*r^2 - 45))/10];
            
    dN( 9,:) = [(t*(128*r*s - 135*t + 45))/10 - (64*s*t*u)/5 ...
                (t*(128*r*s - 135*t + 45))/10 - (64*r*t*u)/5 ...
                (27*t*u)/2 + (t*(128*r*s - 135*t + 45))/10 - (u*(128*r*s - 135*t + 45))/10];
            
    dN(10,:) = [-(t*(270*u - 128*r*s + 128*s*u - 45))/10 ...
                -(t*(270*u - 128*r*s + 128*r*u - 45))/10 ...
                 (t*(128*r*s - 135*u + 45))/10 - (27*t*u)/2 - (u*(128*r*s - 135*u + 45))/10];
            
    dN(11,:) = [ (r*(128*s*t - 135*u + 45))/10 - (27*r*u)/2 - (u*(128*s*t - 135*u + 45))/10 ...
                -(r*(270*u - 128*s*t + 128*t*u - 45))/10 ...
                -(r*(270*u - 128*s*t + 128*s*u - 45))/10];

    dN(12,:) = [(27*r*u)/2 + (r*(128*s*t - 135*r + 45))/10 - (u*(128*s*t - 135*r + 45))/10 ...
                (r*(128*s*t - 135*r + 45))/10 - (64*r*t*u)/5 ...
                (r*(128*s*t - 135*r + 45))/10 - (64*r*s*u)/5];
            
    dN(13,:) = [(t*(270*r - 128*s + 256*r*s + 128*s*t + 128*s^2 - 45))/10 ...
                (64*r*t*(s - u))/5 ...
                (r*(135*r + 128*s*t - 128*s*u - 45))/10];

    dN(14,:) = [(t*(135*t + 128*r*s - 128*s*u - 45))/10 ...
                (64*r*t*(s - u))/5 ...
                (r*(270*t - 128*s + 128*r*s + 256*s*t + 128*s^2 - 45))/10];

    dN(15,:) = [-(s*(270*u - 128*r*t + 128*t*u - 45))/10 ...
                 (s*(128*r*t - 135*u + 45))/10 - (27*s*u)/2 - (u*(128*r*t - 135*u + 45))/10 ...
                -(s*(270*u - 128*r*t + 128*r*u - 45))/10];

    dN(16,:) = [(s*(128*r*t - 135*s + 45))/10 - (64*s*t*u)/5 ...
                (27*s*u)/2 + (s*(128*r*t - 135*s + 45))/10 - (u*(128*r*t - 135*s + 45))/10 ...
                (s*(128*r*t - 135*s + 45))/10 - (64*r*s*u)/5];

    dN(17,:) = [-(s*t*(64*u - 64*r + 135))/5 ...
                 (t*(64*r - 135)*(s - u))/5 ...
                 (s*(64*r - 135)*(t - u))/5];

    dN(18,:) = [ (t*(64*s - 135)*(r - u))/5 ...
                -(r*t*(64*u - 64*s + 135))/5 ...
                 (r*(64*s - 135)*(t - u))/5];

    dN(19,:) = [ (s*(64*t - 135)*(r - u))/5 ...
                 (r*(64*t - 135)*(s - u))/5 ...
                -(r*s*(64*u - 64*t + 135))/5];

    dN(20,:) = [ (s*t*(64*r - 64*u + 135))/5 ...
                 (r*t*(64*s - 64*u + 135))/5 ...
                 (r*s*(64*t - 64*u + 135))/5];

    dN(21,:) = [-256*s*t*(r - u) ...
                -256*r*t*(s - u) ...
                -256*r*s*(t - u)];
end

end % END OF FUNCTION sf_dsf_tet21

% #########################################################################

function [N,dN] = sf_dsf_tet20(r,s,t)
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
%  top be checked
% 
r   = r(:)';
s   = s(:)';
t   = t(:)';
u   = 1-r-s-t;

% SHAPE FUNCTION VALUES
N       = zeros(20,length(r));
N( 1,:) = r.*(3.*r - 2).*((3.*r)./2 - 1/2); % 4 vertex nodes
N( 2,:) = s.*(3.*s - 2).*((3.*s)./2 - 1/2);
N( 3,:) = t.*(3.*t - 2).*((3.*t)./2 - 1/2);
N( 4,:) = u.*(3.*u - 2).*((3.*u)./2 - 1/2);
N( 5,:) = r.*s.*((27.*r)./2 - 9/2); % 12 edge nodes
N( 6,:) = r.*s.*((27.*s)./2 - 9/2);
N( 7,:) = s.*t.*((27.*s)./2 - 9/2);
N( 8,:) = s.*t.*((27.*t)./2 - 9/2);
N( 9,:) = t.*u.*((27.*t)./2 - 9/2);
N(10,:) = t.*u.*((27.*u)./2 - 9/2);
N(11,:) = r.*u.*((27.*u)./2 - 9/2);
N(12,:) = r.*u.*((27.*r)./2 - 9/2);
N(13,:) = r.*t.*((27.*r)./2 - 9/2);
N(14,:) = r.*t.*((27.*t)./2 - 9/2);
N(15,:) = s.*u.*((27.*u)./2 - 9/2);
N(16,:) = s.*u.*((27.*s)./2 - 9/2);
N(17,:) = 27.*s.*t.*u; % 4 face nodes
N(18,:) = 27.*r.*t.*u;
N(19,:) = 27.*r.*s.*u;
N(20,:) = 27.*r.*s.*t;

% DERIVATIVES
if nargout==2
    dN       = zeros(20,3); % derivatives (3 for each node)

    dN( 1,:) = [(27*r^2)/2 - 9*r + 1 ...
                 0 ...
                 0];

    dN( 2,:) = [ 0 ...
                (27*s^2)/2 - 9*s + 1 ...
                 0];

    dN( 3,:) = [ 0 ...
                 0 ...
                (27*t^2)/2 - 9*t + 1];

    dN( 4,:) = [ 25/2 - 27*r*s - 27*r*t - 27*s*t - (27*r^2)/2 - (27*s^2)/2 - (27*t^2)/2 - 18*u ...
                 25/2 - 27*r*s - 27*r*t - 27*s*t - (27*r^2)/2 - (27*s^2)/2 - (27*t^2)/2 - 18*u ...
                 25/2 - 27*r*s - 27*r*t - 27*s*t - (27*r^2)/2 - (27*s^2)/2 - (27*t^2)/2 - 18*u];

    dN( 5,:) = [ (9*s*(6*r - 1))/2 ...
                 (9*r*(3*r - 1))/2 ...
                 0];

    dN( 6,:) = [ (9*s*(3*s - 1))/2 ...
                 (9*r*(6*s - 1))/2 ...
                 0];

    dN( 7,:) = [ 0 ...
                 (9*t*(6*s - 1))/2 ...
                 (9*s*(3*s - 1))/2];

    dN( 8,:) = [ 0 ...
                 (9*t*(3*t - 1))/2 ...
                 (9*s*(6*t - 1))/2];

    dN( 9,:) = [-(9*t*(3*t - 1))/2 ...
                -(9*t*(3*t - 1))/2 ...
                 (9*r)/2 + (9*s)/2 + 36*t - 27*r*t - 27*s*t - (81*t^2)/2 - 9/2];

    dN(10,:) = [-(9*t*(6*u - 1))/2 ...
                -(9*t*(6*u - 1))/2 ...
                 (9*t)/2 + (9*u*(3*u - 1))/2 - 27*t*u];

    dN(11,:) = [ (9*r)/2 + (9*u*(3*u - 1))/2 - 27*r*u ...
                -(9*r*(6*u - 1))/2 ...
                -(9*r*(6*u - 1))/2];

    dN(12,:) = [ 36*r + (9*s)/2 + (9*t)/2 - 27*r*s - 27*r*t - (81*r^2)/2 - 9/2 ...
                -(9*r*(3*r - 1))/2 ...
                -(9*r*(3*r - 1))/2];

    dN(13,:) = [(9*t*(6*r - 1))/2 ...
                 0 ...
                (9*r*(3*r - 1))/2];

    dN(14,:) = [(9*t*(3*t - 1))/2 ...
                 0 ...
                (9*r*(6*t - 1))/2];

    dN(15,:) = [-(9*s*(6*u - 1))/2 ...
                 (9*s)/2 + (9*u*(3*u - 1))/2 - 27*s*u ...
                -(9*s*(6*u - 1))/2];

    dN(16,:) = [-(9*s*(3*s - 1))/2 ...
                 (9*r)/2 + 36*s + (9*t)/2 - 27*r*s - 27*s*t - (81*s^2)/2 - 9/2 ...
                -(9*s*(3*s - 1))/2];

    dN(17,:) = [-27*s*t ...
                -27*t*(s - u) ...
                -27*s*(t - u)];

    dN(18,:) = [-27*t*(r - u) ...
                -27*r*t ...
                -27*r*(t - u)];

    dN(19,:) = [-27*s*(r - u) ...
                -27*r*(s - u) ...
                -27*r*s];

    dN(20,:) = [ 27*s*t ...
                 27*r*t ...
                 27*r*s];
end

end % END OF FUNCTION sf_dsf_tet20

% #########################################################################

% function [N,dN] = sf_dsf_tet15(r,s,t)
% % Find shape functions and their derivatives at given points on the
% % master element for a 15 node tetrahedron
% 
% % Note that shape functions 1:10 sum up to one and 11:15 are added to
% % these. That means values at nodes 11:15 are corrections to the velocity
% % defined by shape functions 1:10!
% 
% % Node notation taken from Hughes' book, p.171;
% % see also Zienkiewicz (4th) Volume 1, p.137
% % local coords of the 15 nodes:
% % Node  r   s   t
% %  1    1   0   0
% %  2    0   1   0
% %  3    0   0   1
% %  4    0   0   0
% %  5   .5  .5   0
% %  6    0  .5  .5
% %  7    0   0  .5
% %  8   .5   0   0
% %  9   .5   0  .5
% % 10    0  .5   0
% % 11    0  1/3 1/3
% % 12   1/3  0  1/3
% % 13   1/3 1/3  0
% % 14   1/3 1/3 1/3
% % 15   .25 .25 .25
% 
% r   = r(:)';
% s   = s(:)';
% t   = t(:)';
% u   = 1-r-s-t;
% 
% % SHAPE FUNCTION VALUES
% N       = zeros(15,length(r));
% N( 1,:) = r.*(2.*r - 1); % 4 vertex nodes
% N( 2,:) = s.*(2.*s - 1);
% N( 3,:) = t.*(2.*t - 1);
% N( 4,:) = u.*(2.*u - 1);
% N( 5,:) =      4.*r.*s; % 6 edge nodes
% N( 6,:) =      4.*s.*t;
% N( 7,:) =      4.*t.*u;
% N( 8,:) =      4.*r.*u;
% N( 9,:) =      4.*r.*t;
% N(10,:) =      4.*s.*u;
% N(11,:) =  27.*s.*t.*u; % 4 face nodes
% N(12,:) =  27.*r.*t.*u;
% N(13,:) =  27.*r.*s.*u;
% N(14,:) =  27.*r.*s.*t;
% N(15,:) = 256.*r.*s.*t.*u; % central buuble
%  
% % DERIVATIVES
% if nargout==2
%     dN       = zeros(15,3); % derivatives (3 for each node)
% 
%     dN( 1,:) = [4*r - 1 ...
%                 0 ...
%                 0];
% 
%     dN( 2,:) = [0 ...
%                 4*s - 1 ...
%                 0];
% 
%     dN( 3,:) = [0 ...
%                 0 ...
%                 4*t - 1];
% 
%     dN( 4,:) = [1 - 4*u ...
%                 1 - 4*u ...
%                 1 - 4*u];
% 
%     dN( 5,:) = [4*s ...
%                 4*r ...
%                 0];
% 
%     dN( 6,:) = [0 ...
%                 4*t ...
%                 4*s];
% 
%     dN( 7,:) = [-4*t ...
%                 -4*t ...
%                 4*u - 4*t];
% 
%     dN( 8,:) = [4*u - 4*r ...
%                 -4*r ...
%                 -4*r];
% 
%     dN( 9,:) = [4*t ...
%                 0 ...
%                 4*r];
% 
%     dN(10,:) = [-4*s ...
%                 4*u - 4*s ...
%                 -4*s];
% 
%     dN(11,:) = [-27*s*t ...
%                 -27*t*(s - u) ...
%                 -27*s*(t - u)];
% 
%     dN(12,:) = [-27*t*(r - u) ...
%                 -27*r*t ...
%                 -27*r*(t - u)];
% 
%     dN(13,:) = [-27*s*(r - u) ...
%                 -27*r*(s - u) ...
%                 -27*r*s];
% 
%     dN(14,:) = [27*s*t ...
%                 27*r*t ...
%                 27*r*s];
% 
%     dN(15,:) = [-256*s*t*(r - u) ...
%                 -256*r*t*(s - u) ...
%                 -256*r*s*(t - u)];
% end
% 
% end % END OF FUNCTION sf_dsf_tet15

% #########################################################################

function [N,dN] = sf_dsf_tet15(r,s,t)
% Find shape functions and their derivatives at given points on the
% master element for a 15 node tetrahedron

% All shape functions sum up to one everywhere inside the tetrahedron.

% Node notation taken from Hughes' book, p.171;
% see also Zienkiewicz (4th) Volume 1, p.137
% local coords of the 15 nodes:
% Node  r   s   t
%  1    1   0   0
%  2    0   1   0
%  3    0   0   1
%  4    0   0   0
%  5   .5  .5   0
%  6    0  .5  .5
%  7    0   0  .5
%  8   .5   0   0
%  9   .5   0  .5
% 10    0  .5   0
% 11    0  1/3 1/3
% 12   1/3  0  1/3
% 13   1/3 1/3  0
% 14   1/3 1/3 1/3
% 15   .25 .25 .25

r   = r(:);
s   = s(:);
t   = t(:);
u   = 1-r-s-t;

% SHAPE FUNCTION VALUES
N = [r.*(2.*r - 1) - 27.*(r.*s.*t + r.*s.*u + r.*t.*u + s.*t.*u)./10 - (128.*r.*s.*t.*u)./7 ... % 4 vertex nodes
     s.*(2.*s - 1) - 27.*(r.*s.*t + r.*s.*u + r.*t.*u + s.*t.*u)./10 - (128.*r.*s.*t.*u)./7 ...
     t.*(2.*t - 1) - 27.*(r.*s.*t + r.*s.*u + r.*t.*u + s.*t.*u)./10 - (128.*r.*s.*t.*u)./7 ...
     u.*(2.*u - 1) - 27.*(r.*s.*t + r.*s.*u + r.*t.*u + s.*t.*u)./10 - (128.*r.*s.*t.*u)./7 ...
     4.*r.*s       - 27.*(r.*s.*t + r.*s.*u + r.*t.*u + s.*t.*u)./10 - (128.*r.*s.*t.*u)./7 ... % 6 edge nodes
     4.*s.*t       - 27.*(r.*s.*t + r.*s.*u + r.*t.*u + s.*t.*u)./10 - (128.*r.*s.*t.*u)./7 ...
     4.*t.*u       - 27.*(r.*s.*t + r.*s.*u + r.*t.*u + s.*t.*u)./10 - (128.*r.*s.*t.*u)./7 ...
     4.*r.*u       - 27.*(r.*s.*t + r.*s.*u + r.*t.*u + s.*t.*u)./10 - (128.*r.*s.*t.*u)./7 ...
     4.*r.*t       - 27.*(r.*s.*t + r.*s.*u + r.*t.*u + s.*t.*u)./10 - (128.*r.*s.*t.*u)./7 ...
     4.*s.*u       - 27.*(r.*s.*t + r.*s.*u + r.*t.*u + s.*t.*u)./10 - (128.*r.*s.*t.*u)./7 ...
     27.*s.*t.*u   - (128.*r.*s.*t.*u)./7 ... % 4 face nodes
     27.*r.*t.*u   - (128.*r.*s.*t.*u)./7 ...
     27.*r.*s.*u   - (128.*r.*s.*t.*u)./7 ...
     27.*r.*s.*t   - (128.*r.*s.*t.*u)./7 ...
     256.*r.*s.*t.*u]'; % central bubble
 
% DERIVATIVES
if nargout==2
    dN       = zeros(15,3); % derivatives (3 for each node)

    dN( 1,:) = [4*r - (27*s)/10 - (27*t)/10 + (27*r*s)/5 + (27*r*t)/5 - (451*s*t)/35 + (128*s*t^2)/7 + (128*s^2*t)/7 + (27*s^2)/10 + (27*t^2)/10 + (256*r*s*t)/7 - 1 ...
                ((s - u)*(189*r + 189*t + 1280*r*t))/70 ...
                ((t - u)*(189*r + 189*s + 1280*r*s))/70];

    dN( 2,:) = [((r - u)*(189*s + 189*t + 1280*s*t))/70 ...
                4*s - (27*r)/10 - (27*t)/10 + (27*r*s)/5 - (451*r*t)/35 + (27*s*t)/5 + (128*r*t^2)/7 + (128*r^2*t)/7 + (27*r^2)/10 + (27*t^2)/10 + (256*r*s*t)/7 - 1 ...
                ((t - u)*(189*r + 189*s + 1280*r*s))/70];

    dN( 3,:) = [((r - u)*(189*s + 189*t + 1280*s*t))/70 ...
                ((s - u)*(189*r + 189*t + 1280*r*t))/70 ...
                4*t - (27*s)/10 - (27*r)/10 - (451*r*s)/35 + (27*r*t)/5 + (27*s*t)/5 + (128*r*s^2)/7 + (128*r^2*s)/7 + (27*r^2)/10 + (27*s^2)/10 + (256*r*s*t)/7 - 1];

    dN( 4,:) = [4*r + (13*s)/10 + (13*t)/10 + (27*r*s)/5 + (27*r*t)/5 - (451*s*t)/35 + (128*s*t^2)/7 + (128*s^2*t)/7 + (27*s^2)/10 + (27*t^2)/10 + (256*r*s*t)/7 - 3 ...
                (13*r)/10 + 4*s + (13*t)/10 + (27*r*s)/5 - (451*r*t)/35 + (27*s*t)/5 + (128*r*t^2)/7 + (128*r^2*t)/7 + (27*r^2)/10 + (27*t^2)/10 + (256*r*s*t)/7 - 3 ...
                (13*r)/10 + (13*s)/10 + 4*t - (451*r*s)/35 + (27*r*t)/5 + (27*s*t)/5 + (128*r*s^2)/7 + (128*r^2*s)/7 + (27*r^2)/10 + (27*s^2)/10 + (256*r*s*t)/7 - 3];

    dN( 5,:) = [(13*s)/10 - (27*t)/10 + (27*r*s)/5 + (27*r*t)/5 - (451*s*t)/35 + (128*s*t^2)/7 + (128*s^2*t)/7 + (27*s^2)/10 + (27*t^2)/10 + (256*r*s*t)/7 ...
                (13*r)/10 - (27*t)/10 + (27*r*s)/5 - (451*r*t)/35 + (27*s*t)/5 + (128*r*t^2)/7 + (128*r^2*t)/7 + (27*r^2)/10 + (27*t^2)/10 + (256*r*s*t)/7 ...
                ((t - u)*(189*r + 189*s + 1280*r*s))/70];

    dN( 6,:) = [((r - u)*(189*s + 189*t + 1280*s*t))/70 ...
                (13*t)/10 - (27*r)/10 + (27*r*s)/5 - (451*r*t)/35 + (27*s*t)/5 + (128*r*t^2)/7 + (128*r^2*t)/7 + (27*r^2)/10 + (27*t^2)/10 + (256*r*s*t)/7 ...
                (13*s)/10 - (27*r)/10 - (451*r*s)/35 + (27*r*t)/5 + (27*s*t)/5 + (128*r*s^2)/7 + (128*r^2*s)/7 + (27*r^2)/10 + (27*s^2)/10 + (256*r*s*t)/7];

    dN( 7,:) = [(27*r*s)/5 - (67*t)/10 - (27*s)/10 + (27*r*t)/5 - (451*s*t)/35 + (128*s*t^2)/7 + (128*s^2*t)/7 + (27*s^2)/10 + (27*t^2)/10 + (256*r*s*t)/7 ...
                (27*r*s)/5 - (67*t)/10 - (27*r)/10 - (451*r*t)/35 + (27*s*t)/5 + (128*r*t^2)/7 + (128*r^2*t)/7 + (27*r^2)/10 + (27*t^2)/10 + (256*r*s*t)/7 ...
                ((t - u)*(189*r + 189*s + 1280*r*s - 280))/70];

    dN( 8,:) = [((r - u)*(189*s + 189*t + 1280*s*t - 280))/70 ...
                (27*r*s)/5 - (27*t)/10 - (67*r)/10 - (451*r*t)/35 + (27*s*t)/5 + (128*r*t^2)/7 + (128*r^2*t)/7 + (27*r^2)/10 + (27*t^2)/10 + (256*r*s*t)/7 ...
                (27*r*t)/5 - (27*s)/10 - (451*r*s)/35 - (67*r)/10 + (27*s*t)/5 + (128*r*s^2)/7 + (128*r^2*s)/7 + (27*r^2)/10 + (27*s^2)/10 + (256*r*s*t)/7];

    dN( 9,:) = [(13*t)/10 - (27*s)/10 + (27*r*s)/5 + (27*r*t)/5 - (451*s*t)/35 + (128*s*t^2)/7 + (128*s^2*t)/7 + (27*s^2)/10 + (27*t^2)/10 + (256*r*s*t)/7 ...
                ((s - u)*(189*r + 189*t + 1280*r*t))/70 ...
                (13*r)/10 - (27*s)/10 - (451*r*s)/35 + (27*r*t)/5 + (27*s*t)/5 + (128*r*s^2)/7 + (128*r^2*s)/7 + (27*r^2)/10 + (27*s^2)/10 + (256*r*s*t)/7];

    dN(10,:) = [(27*r*s)/5 - (27*t)/10 - (67*s)/10 + (27*r*t)/5 - (451*s*t)/35 + (128*s*t^2)/7 + (128*s^2*t)/7 + (27*s^2)/10 + (27*t^2)/10 + (256*r*s*t)/7 ...
                ((s - u)*(189*r + 189*t + 1280*r*t - 280))/70 ...
                (27*r*t)/5 - (67*s)/10 - (451*r*s)/35 - (27*r)/10 + (27*s*t)/5 + (128*r*s^2)/7 + (128*r^2*s)/7 + (27*r^2)/10 + (27*s^2)/10 + (256*r*s*t)/7];

    dN(11,:) = [-(s*t*(128*u - 128*r + 189))/7 ...
                (t*(128*r - 189)*(s - u))/7 ...
                (s*(128*r - 189)*(t - u))/7];

    dN(12,:) = [(t*(128*s - 189)*(r - u))/7 ...
                -(r*t*(128*u - 128*s + 189))/7 ...
                (r*(128*s - 189)*(t - u))/7];

    dN(13,:) = [(s*(128*t - 189)*(r - u))/7 ...
                (r*(128*t - 189)*(s - u))/7 ...
                -(r*s*(128*u - 128*t + 189))/7];

    dN(14,:) = [(s*t*(128*r - 128*u + 189))/7 ...
                (r*t*(128*s - 128*u + 189))/7 ...
                (r*s*(128*t - 128*u + 189))/7];

    dN(15,:) = [-256*s*t*(r - u) ...
                -256*r*t*(s - u) ...
                -256*r*s*(t - u)];
end

end % END OF FUNCTION sf_dsf_tet15

% #########################################################################

function [N,dN] = sf_dsf_tet14(r,s,t)
% Find shape functions and their derivatives at given points on the
% master element for a 14 node tetrahedron

% All shape functions sum up to one everywhere inside the tetrahedron.

% Node notation taken from Hughes' book, p.171;
% see also Zienkiewicz (4th) Volume 1, p.137
% local coords of the 14 nodes:
% Node  r   s   t
%  1    1   0   0
%  2    0   1   0
%  3    0   0   1
%  4    0   0   0
%  5   .5  .5   0
%  6    0  .5  .5
%  7    0   0  .5
%  8   .5   0   0
%  9   .5   0  .5
% 10    0  .5   0
% 11    0  1/3 1/3
% 12   1/3  0  1/3
% 13   1/3 1/3  0
% 14   1/3 1/3 1/3

r   = r(:);
s   = s(:);
t   = t(:);
u   = 1-r-s-t;

% SHAPE FUNCTION VALUES                                                                    
N = [r.*(2.*r - 1) + 27.*(          r.*t.*u + r.*s.*u + r.*s.*t)./9 ... % 4 vertex nodes
     s.*(2.*s - 1) + 27.*(s.*t.*u           + r.*s.*u + r.*s.*t)./9 ...
     t.*(2.*t - 1) + 27.*(s.*t.*u + r.*t.*u 	      + r.*s.*t)./9 ...
     u.*(2.*u - 1) + 27.*(s.*t.*u + r.*t.*u + r.*s.*u          )./9 ...
     4.*r.*s       - 27.*(                    r.*s.*u + r.*s.*t).*4/9 ... % 6 edge nodes
     4.*s.*t       - 27.*(s.*t.*u                     + r.*s.*t).*4/9 ...
     4.*t.*u       - 27.*(s.*t.*u + r.*t.*u                    ).*4/9 ...
     4.*r.*u       - 27.*(          r.*t.*u + r.*s.*u          ).*4/9 ...
     4.*r.*t       - 27.*(          r.*t.*u           + r.*s.*t).*4/9 ...
     4.*s.*u       - 27.*(s.*t.*u           + r.*s.*u          ).*4/9 ...
     27.*s.*t.*u   ... % 4 face nodes
     27.*r.*t.*u   ...
     27.*r.*s.*u   ...
     27.*r.*s.*t]';

% DERIVATIVES
if nargout==2
    dN       = zeros(14,3); % derivatives (3 for each node)

    dN( 1,:) = [4*r - 1 + 3*(s * (u - r) + t * (u - r + s)) ...
                          3*(r * (u - s))                   ...
                          3*(r * (u - t))];
    

    dN( 2,:) = [          3*(s * (u - r))                   ...
                4*s - 1 + 3*(r * (u - s) + t * (u + r - s)) ...
                          3*(s * (u - t))];

    dN( 3,:) = [          3*(t * (u - r))                   ...
                          3*(t * (u - s))                   ...
                4*t - 1 + 3*(r * (u - t) + s * (u + r - t))];

    dN( 4,:) = [4*r + 4*s + 4*t - 3 + 3*(s * (u - r) + t * (u - r - s)) ...
                4*r + 4*s + 4*t - 3 + 3*(r * (u - s) + t * (u - r - s)) ...
                4*r + 4*s + 4*t - 3 + 3*(r * (u - t) + s * (u - r - t))];

    dN( 5,:) = [4*s * (1 - 3*(u + t - r)) ...
                4*r * (1 - 3*(u + t - s)) ...
                0];

    dN( 6,:) = [0                         ...
                4*t * (1 - 3*(u + r - s)) ...
                4*s * (1 - 3*(u + r - t))];

    dN( 7,:) = [-4*t * (1 + 3*(u - r - s)) ...
                -4*t * (1 + 3*(u - r - s)) ...
                 4*(u - t) * (1 - 3*(r + s))];

    dN( 8,:) = [ 4*(u - r) * (1 - 3*(s + t)) ...
                -4*r * (1 + 3*(u - s - t))   ...
                -4*r * (1 + 3*(u - s - t))];

    dN( 9,:) = [4*t * (1 - 3*(u + s - r)) ...
                0                         ...
                4*r * (1 - 3*(u + s - t))];

    dN(10,:) = [-4*s * (1 + 3*(u - r - t))   ...
                 4*(u - s) * (1 - 3*(r + t)) ...
                -4*s * (1 + 3*(u - r - t))];

    dN(11,:) = [-27*s*t ...
                 27*t*(u - s) ...
                 27*s*(u - t)];

    dN(12,:) = [ 27*t*(u - r) ...
                -27*r*t ...
                 27*r*(u - t)];

    dN(13,:) = [ 27*s*(u - r) ...
                 27*r*(u - s) ...
                -27*r*s];

    dN(14,:) = [27*s*t ...
                27*r*t ...
                27*r*s];
end

end % END OF FUNCTION sf_dsf_tet14

% #########################################################################

function [N,dN] = sf_dsf_tet11(r,s,t)
% Find shape functions and their derivatives at given points on the
% master element for a 11 node tetrahedron

% All shape functions sum up to one everywhere inside the tetrahedron.

% Node notation taken from Hughes' book, p.171;
% see also Zienkiewicz (4th) Volume 1, p.137
% local coords of the 11 nodes:
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
% 11  .25 .25 .25

r   = r(:)';
s   = s(:)';
t   = t(:)';
u   = 1-r-s-t;

% SHAPE FUNCTION VALUES
N       = zeros(11,length(r));
N( 1,:) = r.*(2.*r - 1) - (128.*r.*s.*t.*u)./5; % 4 vertex nodes
N( 2,:) = s.*(2.*s - 1) - (128.*r.*s.*t.*u)./5;
N( 3,:) = t.*(2.*t - 1) - (128.*r.*s.*t.*u)./5;
N( 4,:) = u.*(2.*u - 1) - (128.*r.*s.*t.*u)./5;
N( 5,:) = 4.*r.*s - (128.*r.*s.*t.*u)./5; % 6 edge nodes
N( 6,:) = 4.*s.*t - (128.*r.*s.*t.*u)./5;
N( 7,:) = 4.*t.*u - (128.*r.*s.*t.*u)./5;
N( 8,:) = 4.*r.*u - (128.*r.*s.*t.*u)./5;
N( 9,:) = 4.*r.*t - (128.*r.*s.*t.*u)./5;
N(10,:) = 4.*s.*u - (128.*r.*s.*t.*u)./5;
N(11,:) =            256.*r.*s.*t.*u; % central buuble

% DERIVATIVES
if nargout==2
    dN       = zeros(11,3); % derivatives (3 for each node)

    dN( 1,:) = [ 4*r + (128*r*s*t)/5 - (128*s*t*u)/5 - 1 ...
                                     (128*r*t*(s - u))/5 ...
                                     (128*r*s*(t - u))/5];
                 
    dN( 2,:) = [                     (128*s*t*(r - u))/5 ...
                 4*s + (128*r*s*t)/5 - (128*r*t*u)/5 - 1 ...
                                     (128*r*s*(t - u))/5];
                 
    dN( 3,:) = [                     (128*s*t*(r - u))/5 ...
                                     (128*r*t*(s - u))/5 ...
                 4*t + (128*r*s*t)/5 - (128*r*s*u)/5 - 1];

    dN( 4,:) = [ (128*r*s*t)/5 - 4*u - (128*s*t*u)/5 + 1 ...
                 (128*r*s*t)/5 - 4*u - (128*r*t*u)/5 + 1 ...
                 (128*r*s*t)/5 - 4*u - (128*r*s*u)/5 + 1];

    dN( 5,:) = [ (4*s*(64*r*t - 32*t + 32*s*t + 32*t^2 + 5))/5 ...
                 (4*r*(32*r*t - 32*t + 64*s*t + 32*t^2 + 5))/5 ...
                                           (128*r*s*(t - u))/5];
                       
    dN( 6,:) = [                           (128*s*t*(r - u))/5 ...
                 (4*t*(64*r*s - 32*r + 32*r*t + 32*r^2 + 5))/5 ...
                 (4*s*(32*r*s - 32*r + 64*r*t + 32*r^2 + 5))/5];

    dN( 7,:) = [-(4*t*(32*s*u - 32*r*s + 5))/5 ...
                -(4*t*(32*r*u - 32*r*s + 5))/5 ...
                 (4*(32*r*s - 5)*(t - u))/5];
 
    dN( 8,:) = [ (4*(r - u)*(32*s*t - 5))/5 ...
                -(4*r*(32*t*u - 32*s*t + 5))/5 ...
                -(4*r*(32*s*u - 32*s*t + 5))/5];

    dN( 9,:) = [ (4*t*(64*r*s - 32*s + 32*s*t + 32*s^2 + 5))/5 ...
                                           (128*r*t*(s - u))/5 ...
                 (4*r*(32*r*s - 32*s + 64*s*t + 32*s^2 + 5))/5];

    dN(10,:) = [-(4*s*(32*t*u - 32*r*t + 5))/5 ...
                 (4*(32*r*t - 5)*(s - u))/5 ...
                -(4*s*(32*r*u - 32*r*t + 5))/5];

    dN(11,:) = [-256*s*t*(r - u) ...
                -256*r*t*(s - u) ...
                -256*r*s*(t - u)];
end

end % END OF FUNCTION sf_dsf_tet11

% #########################################################################

function [N,dN] = sf_dsf_tet10(r,s,t)
% Find shape functions and their derivatives at given points on the
% master element for a 10 node tetrahedron

% All shape functions sum up to one everywhere inside the tetrahedron.

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

r   = r(:)';
s   = s(:)';
t   = t(:)';
u   = 1-r-s-t;

% SHAPE FUNCTION VALUES
N       = zeros(10,length(r));
N( 1,:) = r.*(2.*r-1); % 4 vertex nodes
N( 2,:) = s.*(2.*s-1);
N( 3,:) = t.*(2.*t-1);
N( 4,:) = u.*(2.*u-1);
N( 5,:) = 4.*r.*s; % 6 edge nodes
N( 6,:) = 4.*s.*t;
N( 7,:) = 4.*t.*u;
N( 8,:) = 4.*r.*u;
N( 9,:) = 4.*r.*t;
N(10,:) = 4.*s.*u;

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

    dN(4,1) = 1-4*u; % dN/dr
    dN(4,2) = 1-4*u; % dN/ds
    dN(4,3) = 1-4*u; % dN/dt

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

% All shape functions sum up to one everywhere inside the tetrahedron.

% Node notation taken from Hughes' book, p.170
% local coords of the 4 nodes:
% Node  r  s  t
%  1    1  0  0
%  2    0  1  0
%  3    0  0  1
%  4    0  0  0

r  = r(:)';
s  = s(:)';
t  = t(:)';
u  = 1-r-s-t;

% SHAPE FUNCTION VALUES
N  = [r; s; t; u];

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