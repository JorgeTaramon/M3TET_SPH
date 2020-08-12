function len_el = calc_tetra_length(GCOORD,EL2NOD)
%
% Purpose: Calculate characteristic length for each tetrahedra in the mesh
%
% Input:
%   GCOORD   : [matrix] : coordinates of all nodes in mesh (3 x nnod)
%   EL2NOD   : [matrix] : connectivity matrix (nnodel x nel)
%
% Output:
%   len_el   : [vector] : length of each element
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

nel       = size(EL2NOD,2);

x_el      = reshape( GCOORD(1,EL2NOD(1:4,:)), 4, nel );
dx_6edges = [ ( x_el(1,:) - x_el(2,:) ).^2 ;
              ( x_el(1,:) - x_el(3,:) ).^2 ;
              ( x_el(1,:) - x_el(4,:) ).^2 ;
              ( x_el(2,:) - x_el(3,:) ).^2 ;
              ( x_el(2,:) - x_el(4,:) ).^2 ;
              ( x_el(3,:) - x_el(4,:) ).^2 ];
clear x_el

y_el 	  = reshape( GCOORD(2,EL2NOD(1:4,:)), 4, nel );
dy_6edges = [ ( y_el(1,:) - y_el(2,:) ).^2 ;
              ( y_el(1,:) - y_el(3,:) ).^2 ;
              ( y_el(1,:) - y_el(4,:) ).^2 ;
              ( y_el(2,:) - y_el(3,:) ).^2 ;
              ( y_el(2,:) - y_el(4,:) ).^2 ;
              ( y_el(3,:) - y_el(4,:) ).^2 ];
clear y_el

z_el      = reshape( GCOORD(3,EL2NOD(1:4,:)), 4, nel );
dz_6edges = [ ( z_el(1,:) - z_el(2,:) ).^2 ;
              ( z_el(1,:) - z_el(3,:) ).^2 ;
              ( z_el(1,:) - z_el(4,:) ).^2 ;
              ( z_el(2,:) - z_el(3,:) ).^2 ;
              ( z_el(2,:) - z_el(4,:) ).^2 ;
              ( z_el(3,:) - z_el(4,:) ).^2 ];
clear z_el

length_6edges = sqrt( dx_6edges + dy_6edges + dz_6edges );
len_el        = min(length_6edges)';

end % END OF SUBFUNCTION tetra_length