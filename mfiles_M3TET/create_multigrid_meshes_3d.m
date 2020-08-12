function [GCOORD,EL2NOD_mg,PointID,PhaseID,Ic2f] = create_multigrid_meshes_3d...
             (GCOORD_c,EL2NOD_c,PointID_c,PhaseID_c,DB_indices,nmg)
% Usage: [GCOORD,EL2NOD_mg,PointID,PhaseID,Ic2f] = create_multigrid_meshes_3d...
%            (GCOORD_c,EL2NOD_c,PointID_c,PhaseID_c,DB_indices,nmg)
%
% Purpose: Creates "nmg"-1 nested multigrid meshes by recursively refining
%          the input mesh (3D tetrahedra mesh, either 4-node linear or
%          10-node quadratic order elements)
%
% Input:
%   GCOORD_c   : [matrix] : coordinates of all nodes in coarse mesh
%   EL2NOD_c   : [matrix] : connectivity matrix of coarse mesh (nnodel x nel)
%   PointID_c  : [vector] : domain boundary index for each node
%   PhaseID_c  : [vector] : element material index for each element
%   DB_indices : [matrix] : domain boundary indices for each face of the domain
%   nmg        : [scalar] : number of multigrid levels
%
% Output:
%   GCOORD     : [matrix] : coordinates of all nodes in finest mesh
%   EL2NOD_mg  : [cell]   : connectivity matrices for each multigrid level
%   PointID    : [vector] : domain boundary index for each node in finest mesh
%   PhaseID    : [vector] : element material index for each element in finest mesh
%   Ic2f       : [cell]   : transfer operators between meshes
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Mar 2013
% JH Feb 2016 : restructured, cleaned up
%

% Recursively refine the mesh by splitting each tetrahedron into 8 elements
% Nodal coordinates, connectivity matrix and boundary node information are
% recalculated during the refinement. For each refinement, an interpolation
% matrix "Ic2f" is calculated that will be used to transfer data between
% the levels.
EL2NOD_mg      = cell(1,nmg);
EL2NOD_mg{nmg} = uint32(EL2NOD_c);
Ic2f           = cell(1,nmg-1);
ndof           = 3; % number of degrees of freedom per node (velocity --> 3)

for img=nmg:-1:2
    % CALCULATE REFINED MESH
    [GCOORD,EL2NOD,PointID,PhaseID] = tetmesh_refine...
        (GCOORD_c,EL2NOD_c,PointID_c,PhaseID_c,DB_indices);
    
    % STORE CONNECTIVITY MATRIX
    EL2NOD_mg{img-1}        = uint32(EL2NOD);
    
    % CALCULATE INTERMESH-TRANSFER OPERATOR (INTERPOLATION MATRIX)
    Ic2f{img-1} = intermesh_transfer_matrix...
        (GCOORD,EL2NOD,GCOORD_c,EL2NOD_c,ndof);

    % SAVE CURRENT MESH AS 'COARSE'
    if img>2
        GCOORD_c  = GCOORD;
        EL2NOD_c  = EL2NOD;
        PointID_c = PointID;
        PhaseID_c = PhaseID;
    end
end

end % END OF FUNCTION create_multigrid_meshes_3d