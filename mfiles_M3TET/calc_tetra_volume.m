function vol_el = calc_tetra_volume(GCOORD,EL2NOD)

% Purpose: Calculate volume of each tetrahedron in the mesh
%
% Input:
%   GCOORD   : [matrix] : coordinates of all nodes in mesh (3 x nnod)
%   EL2NOD   : [matrix] : connectivity matrix (nnodel x nel)
%
% Output:
%   vol_el   : [vector] : volume of each element
%
% Part of 3D convection code M3TET
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de
% For numerical methods see online Ph.D. thesis
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Mar 2011
% JH Feb 2016 : cleaned up, cleared variables to save memory
%

nel    = size(EL2NOD,2);
nnodel = 4;

% Derivatives of the 4 linear tetrahedra shape functions
dNds = [ 1  0  0   % dN1/dr dN1/ds dN1/dt
         0  1  0   % dN2/dr dN2/ds dN2/dt
         0  0  1   % dN3/dr dN3/ds dN3/dt
        -1 -1 -1]; % dN4/dr dN4/ds dN4/dt

% Components of Jacobi matrices of all elements in block
x_el = reshape( GCOORD(1,EL2NOD(1:4,:)), nnodel, nel );
Jx   = x_el'*dNds; clear x_el
y_el = reshape( GCOORD(2,EL2NOD(1:4,:)), nnodel, nel );
Jy   = y_el'*dNds; clear y_el
z_el = reshape( GCOORD(3,EL2NOD(1:4,:)), nnodel, nel );
Jz   = z_el'*dNds; clear z_el

% Determinats of Jaciobi matrices ("Js")
detJ   =     Jx(:,1).*Jy(:,2).*Jz(:,3) ...
           + Jx(:,2).*Jy(:,3).*Jz(:,1) ...
           + Jx(:,3).*Jy(:,1).*Jz(:,2) ...
           - Jx(:,3).*Jy(:,2).*Jz(:,1) ...
           - Jx(:,1).*Jy(:,3).*Jz(:,2) ...
           - Jx(:,2).*Jy(:,1).*Jz(:,3);
vol_el = detJ ./ 6;

end % END OF FUNCTION calc_tetra_volume