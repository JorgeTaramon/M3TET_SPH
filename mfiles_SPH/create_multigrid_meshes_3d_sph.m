function [GCOORD,GCOORD_SPH,EL2NOD_mg,PointID,PhaseID,Ic2f, ...
    els_out_cone_no_cross_2pi,els_out_cone_cross_2pi,       ...
    els_in_cone_no_iso,els_in_cone_iso,                     ...
    GCOORD_face_nodes,EL2NOD14,GCOORD_cubic,EL2NOD_cubic] = ...
    create_multigrid_meshes_3d_sph(GCOORD_c,GCOORD_SPH_c,   ...
        EL2NOD_c,PointID_c,PhaseID_c,DB_indices,nmg,r_cmb,r_surf,SETTINGS)
% Usage: [GCOORD,GCOORD_SPH,EL2NOD_mg,PointID,PhaseID,Ic2f, ...
%   els_out_cone_no_cross_2pi,els_out_cone_cross_2pi,       ...
%   els_in_cone_no_iso,els_in_cone_iso,                     ...
%   GCOORD_face_nodes,EL2NOD14,GCOORD_cubic,EL2NOD_cubic] = ...
%   create_multigrid_meshes_3d_sph(GCOORD_c,GCOORD_SPH_c,   ...
%       EL2NOD_c,PointID_c,PhaseID_c,DB_indices,nmg,r_cmb,r_surf,SETTINGS)
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

GCOORD_straight_c = GCOORD_c;
nodE = 5:10;
nodV = [1 2; 2 3; 3 4; 4 1; 1 3; 2 4];
for i=1:length(nodE)
    GCOORD_straight_c(1,EL2NOD_c(nodE(i),:)) = ...
        0.5*sum(reshape(GCOORD_c(1,EL2NOD_c(nodV(i,:),:)),2,[]),1);
    GCOORD_straight_c(2,EL2NOD_c(nodE(i),:)) = ...
        0.5*sum(reshape(GCOORD_c(2,EL2NOD_c(nodV(i,:),:)),2,[]),1);
    GCOORD_straight_c(3,EL2NOD_c(nodE(i),:)) = ...
        0.5*sum(reshape(GCOORD_c(3,EL2NOD_c(nodV(i,:),:)),2,[]),1);
end

for img=nmg:-1:2
    % CALCULATE REFINED MESH
    [GCOORD,GCOORD_SPH,EL2NOD,PointID,PhaseID,              ...
    els_out_cone_no_cross_2pi,els_out_cone_cross_2pi,       ...
    els_in_cone_no_iso,els_in_cone_iso,                     ...
    GCOORD_face_nodes,EL2NOD14,GCOORD_cubic,EL2NOD_cubic] = ...
    tetmesh_refine_sph(GCOORD_c,GCOORD_SPH_c,EL2NOD_c,      ...
         PointID_c,PhaseID_c,DB_indices,r_cmb,r_surf,SETTINGS);
    
    % STORE CONNECTIVITY MATRIX
    EL2NOD_mg{img-1}        = uint32(EL2NOD);

    % CALCULATE INTERMESH-TRANSFER OPERATOR (INTERPOLATION MATRIX)
    GCOORD_straight = tetmesh_refine...
        (GCOORD_straight_c,EL2NOD_c,PointID_c,PhaseID_c,DB_indices);
    Ic2f{img-1} = intermesh_transfer_matrix...
        (GCOORD_straight,EL2NOD,GCOORD_straight_c,EL2NOD_c,ndof);

    % SAVE CURRENT MESH AS 'COARSE'
    if img>2
        GCOORD_straight_c = GCOORD_straight;
        GCOORD_c     = GCOORD;
        GCOORD_SPH_c = GCOORD_SPH;
        EL2NOD_c     = EL2NOD;
        PointID_c    = PointID;
        PhaseID_c    = PhaseID;
    end
end

end % END OF FUNCTION create_multigrid_meshes_3d