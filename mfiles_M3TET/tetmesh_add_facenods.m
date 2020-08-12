function [GCOORD,EL2NOD,FACE2NOD] = tetmesh_add_facenods(GCOORD,EL2NOD)

nnodel = size(EL2NOD,1);
nel    = size(EL2NOD,2);
nnod0  = max(EL2NOD(:));

faces = [2 3 4 ... % face 1
         1 3 4 ... % face 2
         1 2 4 ... % face 3
         1 2 3];   % face 4
     
FACE2NOD = reshape( EL2NOD(faces',:),3,[] )';

% Find the faces that are shared by neighboring elements and return a 
% unique list of faces.
[FACE2NOD,~,ib] = unique_keep_order(FACE2NOD);
nface           = size(FACE2NOD,1);

% element to bar connectivity after doubles have been merged (removed)
EL2FACE = reshape(uint32(ib),4,nel);

% write face node connectivity
EL2NOD(nnodel+1:nnodel+4,:) = nnod0 + EL2FACE;

x_facenods = sum(reshape(GCOORD(1,FACE2NOD'),3,nface))./3;
y_facenods = sum(reshape(GCOORD(2,FACE2NOD'),3,nface))./3;
z_facenods = sum(reshape(GCOORD(3,FACE2NOD'),3,nface))./3;
GCOORD     = [GCOORD [x_facenods; y_facenods; z_facenods]];

end