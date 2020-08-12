function [Nip,dNip,pt,wt] = integration_setup_tet(nnodel,nip)
% Usage: [Nip,dNip,pt,wt] = integration_setup_tet(nnodel,nip)
% 
% Purpose: Defines integration points and weights, returns values of shape
%          functions and their derivatives at integration points (IPs)
%
% Input:
%   nnodel : [scalar]   : number of nodes per tetrahedral element
%   nip    : [scalar]   : number of integration points
%
% Output:
%   Nip  : [matrix]    : shape fct values at IPs (nip x nnodel)
%   dNip : [tensor]    : shape fct derivatives at IPs (3 x nnodel x nip)
%   pt   : [matrix]    : local coordinates of IPs (3 x nip)
%   wt   : [rowvector] : integration weights at each IP
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Dec 2012
% JH March 2013
%

switch nip
    case 1
        % ONE POINT INTEGRATION
        % Taken from Zienkiewicz Vol 1 (4th) edition p.177
        % (also in Hughes p.174)
        nip   = 1;
        pt    = zeros(3,nip); % one point defined by the 3 local coords r, s, t
        pt(1) = 1/4; % r
        pt(2) = 1/4; % s
        pt(3) = 1/4; % t
        wt    = 1;   % weight for integration
        wt    = wt ./ 6; % volume of tetrahedron is 1/6 of encompassing cube
                         % V=0.5*A*h, A=area of base, h=height
    case 4
        % FOUR POINT INTEGRATION
        % Taken from Zienkiewicz Vol 1 (4th) edition p.177 (also in Hughes p.174)
        nip = 4;
        a   = (5+3*sqrt(5))/20; %0.58541020; %
        b   = (5-sqrt(5))/20; %0.13819660; %
        pt      = zeros(3,nip); % 4 points defined by the 3 local coords r, s, t
        pt(1,:) = [ a   b   b   b ]; % r at all points
        pt(2,:) = [ b   a   b   b ]; % s at all points
        pt(3,:) = [ b   b   a   b ]; % t at all points
        wt      = [1/4 1/4 1/4 1/4]; % weights for integration
        wt      = wt ./ 6; % volume of tetrahedron is 1/6 of encompassing cube
                           % V=0.5*A*h, A=area of base, h=height
    case 5
        % FIVE POINT INTEGRATION
        % Taken from Zienkiewicz Vol 1 (4th) edition p.177 (also in Hughes p.174)
        nip     = 5;
        pt      = zeros(3,nip); % 5 points defined by the 3 local coords r, s, t
        pt(1,:) = [ 1/4  1/2  1/6  1/6  1/6]; % r at all points
        pt(2,:) = [ 1/4  1/6  1/2  1/6  1/6]; % s at all points
        pt(3,:) = [ 1/4  1/6  1/6  1/2  1/6]; % t at all points
        wt      = [-4/5 9/20 9/20 9/20 9/20]; % weights for integration
        wt      = wt ./ 6; % volume of tetrahedron is 1/6 of encompassing cube
                           % V=0.5*A*h, A=area of base, h=height
    otherwise
        error('Number of integration points (nip) must be 4 or 10.');
end

[Nip,dNip] = shapefct_tet(nnodel,pt(1,:),pt(2,:),pt(3,:));

end % END OF FUNCTION integration_setup_tet