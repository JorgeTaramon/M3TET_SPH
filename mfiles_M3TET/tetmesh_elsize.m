function [length,volume] = tetmesh_elsize(GCOORD,EL2NOD)

% Purpose: Calculate length and volume of each tetrahedron in the mesh
%
% Input:
%   GCOORD     :: nodal coordinates
%   EL2NOD :: connectivity matrix
% Output:
%   ellength :: length-scale for each element
%
% Part of 3D convection code M3TET
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de
% For numerical methods see online Ph.D. thesis
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Mar 2011
% JH Feb 2015
%

nel  = size(EL2NOD,2);
x_el = reshape( GCOORD(1,EL2NOD(1:4,:)), 4, nel );
y_el = reshape( GCOORD(2,EL2NOD(1:4,:)), 4, nel );
z_el = reshape( GCOORD(3,EL2NOD(1:4,:)), 4, nel );

dx_6edges = [ ( x_el(1,:) - x_el(2,:) ).^2 ;
              ( x_el(1,:) - x_el(3,:) ).^2 ;
              ( x_el(1,:) - x_el(4,:) ).^2 ;
              ( x_el(2,:) - x_el(3,:) ).^2 ;
              ( x_el(2,:) - x_el(4,:) ).^2 ;
              ( x_el(3,:) - x_el(4,:) ).^2 ];
dy_6edges = [ ( y_el(1,:) - y_el(2,:) ).^2 ;
              ( y_el(1,:) - y_el(3,:) ).^2 ;
              ( y_el(1,:) - y_el(4,:) ).^2 ;
              ( y_el(2,:) - y_el(3,:) ).^2 ;
              ( y_el(2,:) - y_el(4,:) ).^2 ;
              ( y_el(3,:) - y_el(4,:) ).^2 ];
dz_6edges = [ ( z_el(1,:) - z_el(2,:) ).^2 ;
              ( z_el(1,:) - z_el(3,:) ).^2 ;
              ( z_el(1,:) - z_el(4,:) ).^2 ;
              ( z_el(2,:) - z_el(3,:) ).^2 ;
              ( z_el(2,:) - z_el(4,:) ).^2 ;
              ( z_el(3,:) - z_el(4,:) ).^2 ];
length_6edges = sqrt( dx_6edges + dy_6edges + dz_6edges );
length        = min(length_6edges)';

% Derivatives of the 4 linear tetrahedra shape functions
[~,dNds] = sf_dsf_tet410([1/3 1/3 1/3]',4,'cell');
dNds     = dNds{1};

% Components of Jacobi matrices of all elements in block
Jx = x_el'*dNds;
Jy = y_el'*dNds;
Jz = z_el'*dNds;

% Determinats of Jaciobi matrices ("Js")
detJ   = Jx(:,1).*Jy(:,2).*Jz(:,3) ...
       + Jx(:,2).*Jy(:,3).*Jz(:,1) ...
       + Jx(:,3).*Jy(:,1).*Jz(:,2) ...
       - Jx(:,3).*Jy(:,2).*Jz(:,1) ...
       - Jx(:,1).*Jy(:,3).*Jz(:,2) ...
       - Jx(:,2).*Jy(:,1).*Jz(:,3);
volume = detJ ./ 6;

end