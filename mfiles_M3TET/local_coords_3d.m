function [lcoord] = local_coords_3d(gX,EL2NOD,els,gX_PT)
% Usage: 
%
% Purpose:
%
% Input:
%
% Output:
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH March 2011
% JH March 2013
%

ns = length(els); % Number of points
x  = reshape(gX(1,EL2NOD(1:4,els)),4,ns);
y  = reshape(gX(2,EL2NOD(1:4,els)),4,ns);
z  = reshape(gX(3,EL2NOD(1:4,els)),4,ns);

xp = gX_PT(1,:);
yp = gX_PT(2,:);
zp = gX_PT(3,:);

lcoord = zeros(3,length(els));

lcoord(1,:) =   (  (x(2,:)-xp)    .*(y(3,:)-yp)    .*(z(4,:)-zp)     + (x(3,:)-xp)    .*(y(4,:)-yp)    .*(z(2,:)-zp)     + (x(4,:)-xp)    .*(y(2,:)-yp)    .*(z(3,:)-zp)     ...
                 - (z(2,:)-zp)    .*(y(3,:)-yp)    .*(x(4,:)-xp)     - (z(3,:)-zp)    .*(y(4,:)-yp)    .*(x(2,:)-xp)     - (z(4,:)-zp)    .*(y(2,:)-yp)    .*(x(3,:)-xp) )   ...
             ./ (  (x(2,:)-x(1,:)).*(y(3,:)-y(1,:)).*(z(4,:)-z(1,:)) + (x(3,:)-x(1,:)).*(y(4,:)-y(1,:)).*(z(2,:)-z(1,:)) + (x(4,:)-x(1,:)).*(y(2,:)-y(1,:)).*(z(3,:)-z(1,:))     ...
                 - (z(2,:)-z(1,:)).*(y(3,:)-y(1,:)).*(x(4,:)-x(1,:)) - (z(3,:)-z(1,:)).*(y(4,:)-y(1,:)).*(x(2,:)-x(1,:)) - (z(4,:)-z(1,:)).*(y(2,:)-y(1,:)).*(x(3,:)-x(1,:)) );

lcoord(2,:) =   (  (x(1,:)-xp)    .*(y(3,:)-yp)    .*(z(4,:)-zp)     + (x(3,:)-xp)    .*(y(4,:)-yp)    .*(z(1,:)-zp)     + (x(4,:)-xp)    .*(y(1,:)-yp)    .*(z(3,:)-zp)     ...
                 - (z(1,:)-zp)    .*(y(3,:)-yp)    .*(x(4,:)-xp)     - (z(3,:)-zp)    .*(y(4,:)-yp)    .*(x(1,:)-xp)     - (z(4,:)-zp)    .*(y(1,:)-yp)    .*(x(3,:)-xp) )   ...
             ./ (  (x(1,:)-x(2,:)).*(y(3,:)-y(2,:)).*(z(4,:)-z(2,:)) + (x(3,:)-x(2,:)).*(y(4,:)-y(2,:)).*(z(1,:)-z(2,:)) + (x(4,:)-x(2,:)).*(y(1,:)-y(2,:)).*(z(3,:)-z(2,:))     ...
                 - (z(1,:)-z(2,:)).*(y(3,:)-y(2,:)).*(x(4,:)-x(2,:)) - (z(3,:)-z(2,:)).*(y(4,:)-y(2,:)).*(x(1,:)-x(2,:)) - (z(4,:)-z(2,:)).*(y(1,:)-y(2,:)).*(x(3,:)-x(2,:)) );
             
lcoord(3,:) =   (  (x(1,:)-xp)    .*(y(2,:)-yp)    .*(z(4,:)-zp)     + (x(2,:)-xp)    .*(y(4,:)-yp)    .*(z(1,:)-zp)     + (x(4,:)-xp)    .*(y(1,:)-yp)    .*(z(2,:)-zp)     ...
                 - (z(1,:)-zp)    .*(y(2,:)-yp)    .*(x(4,:)-xp)     - (z(2,:)-zp)    .*(y(4,:)-yp)    .*(x(1,:)-xp)     - (z(4,:)-zp)    .*(y(1,:)-yp)    .*(x(2,:)-xp) )   ...
             ./ (  (x(1,:)-x(3,:)).*(y(2,:)-y(3,:)).*(z(4,:)-z(3,:)) + (x(2,:)-x(3,:)).*(y(4,:)-y(3,:)).*(z(1,:)-z(3,:)) + (x(4,:)-x(3,:)).*(y(1,:)-y(3,:)).*(z(2,:)-z(3,:))     ...
                 - (z(1,:)-z(3,:)).*(y(2,:)-y(3,:)).*(x(4,:)-x(3,:)) - (z(2,:)-z(3,:)).*(y(4,:)-y(3,:)).*(x(1,:)-x(3,:)) - (z(4,:)-z(3,:)).*(y(1,:)-y(3,:)).*(x(2,:)-x(3,:)) );

end % END OF FUNCTION local_coords_3d