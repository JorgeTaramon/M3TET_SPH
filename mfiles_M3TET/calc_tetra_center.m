function elc = calc_tetra_center(GCOORD,EL2NOD)
%
% Purpose: Calculate center coordinates of each tetrahedron in the mesh
% Input:
%   GCOORD :: nodal coordinates
%   EL2NOD :: connectivity matrix
% Output:
%   elc    :: center coordinates of each element
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH March 2011
%

nel      = size(EL2NOD,2);
elc      = zeros(3,nel);
elc(1,:) = 0.25*sum(reshape(GCOORD(1,EL2NOD)',4,nel));
elc(2,:) = 0.25*sum(reshape(GCOORD(2,EL2NOD)',4,nel));
elc(3,:) = 0.25*sum(reshape(GCOORD(3,EL2NOD)',4,nel));

end % END OF FUNCTION calc_tetra_center