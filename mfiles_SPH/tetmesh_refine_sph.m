function [GCOORD,GCOORD_SPH,EL2NOD,PointID,PhaseID,         ...
    els_out_cone_no_cross_2pi,els_out_cone_cross_2pi,       ...
    els_in_cone_no_iso,els_in_cone_iso,                     ...
    GCOORD_face_nodes,EL2NOD14,GCOORD_cubic,EL2NOD_cubic] = ...
    tetmesh_refine_sph(GCOORD_c,GCOORD_SPH_c,EL2NOD_c,      ...
        PointID_c,PhaseID_c,DB_indices,r_cmb,r_surf,SETTINGS)
% Usage: [GCOORD,GCOORD_SPH,EL2NOD,PointID,PhaseID,         ...
%   els_out_cone_no_cross_2pi,els_out_cone_cross_2pi,       ...
%   els_in_cone_no_iso,els_in_cone_iso,                     ...
%   GCOORD_face_nodes,EL2NOD14,GCOORD_cubic,EL2NOD_cubic] = ...
%   tetmesh_refine_sph(GCOORD_c,GCOORD_SPH_c,EL2NOD_c,      ...
%       PointID_c,PhaseID_c,DB_indices,r_cmb,r_surf,SETTINGS)
%
% Purpose: Refines the input mesh (3D tetrahedra mesh, either 4-node linear or
%          10-node quadratic order elements) by splitting each tetrahedron
%          into 8 elements.
%
% Input:
%   GCOORD_c   : [matrix] : coordinates of all nodes in mesh (3 x nnod)
%   EL2NOD_c   : [matrix] : connectivity matrix of mesh (nnodel x nel)
%   PointID_c  : [vector] : domain boundary index for each node (1 x nnod)
%   PhaseID_c  : [vector] : element material index for each element (1 x nel)
%   DB_indices : [matrix] : domain boundary indices for each face of the domain
%
% Output:
%   GCOORD  : [matrix] : coordinates of all nodes in refined mesh (3 x ?)
%   EL2NOD  : [matrix] : connectivity matrix of refined mesh (nnodel x 8*nel)
%   PointID : [vector] : domain boundary index for each node in refined mesh
%   PhaseID : [vector] : element material index for each element in refined mesh
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Mar 2013
% JH Feb 2016 : cleaned up, made this a seperate m-file
%

FigNo  = 0;  % set to zero to NOT show mesh figures

if FigNo
    sfigure(FigNo);clf;
    GCOORD_c = spherical2cartesian(GCOORD_SPH_c);
    tetramesh(EL2NOD_c(1:4,:)',GCOORD_c',...
              'EdgeColor','b',...
              'FaceAlpha',0.1,'FaceColor','k',...
              'LineStyle','-','LineWidth',1);
end

% Check if the mesh has linear (4-node) or quadratic (10-node) tetrahedral
% elements. If quadratic, split each element into 8 linear sub-elements. 
% =========================================================================
nnodel = size(EL2NOD_c,1);
if nnodel==10
    [EL2NOD_c,PhaseID_c] = tetmesh_p2_to_p1(GCOORD_c,EL2NOD_c,PhaseID_c);
end

% =======================================================================================
% CREATE A 10 NODE MESH WITH CURVED EDGE ELEMENTS
% =======================================================================================
[GCOORD,GCOORD_SPH,EL2NOD,PointID,els_out_cone_no_cross_2pi,     ...
    els_out_cone_cross_2pi,els_in_cone_no_iso,els_in_cone_iso] = ...
    tetmesh_p1_to_p2_sph(GCOORD_c,GCOORD_SPH_c,EL2NOD_c,         ...
                         PointID_c,DB_indices,r_cmb,r_surf,SETTINGS);
if ~isempty(els_in_cone_iso)
    % CREATE QUASI-CUBIC ISOPARAMETRIC ELEMENTS (els crossing the cone)
    [GCOORD_face_nodes,EL2NOD14] = face_nodes_els_in_cone_iso ...
        (GCOORD,GCOORD_SPH_c,GCOORD_SPH,EL2NOD,els_in_cone_iso,SETTINGS);
    % CREATE CUBIC ISOPARAMETRIC ELEMENTS (els crossing the cone)
    [GCOORD_cubic,EL2NOD_cubic] = cubic_els_in_cone_iso ...
        (GCOORD_c,GCOORD,GCOORD_SPH_c,EL2NOD_c,els_in_cone_iso,SETTINGS);
else
    GCOORD_face_nodes = [];
    EL2NOD14          = [];
    GCOORD_cubic      = [];
    EL2NOD_cubic      = [];
end

if nnodel==4
    % Re-connect elements and nodes by creating a linear (4-node) element
    % connectivity matrix (only if input was a linear connectivity)
    % =========================================================================
    [EL2NOD,PhaseID] = tetmesh_p2_to_p1(GCOORD_c,EL2NOD_c,PhaseID_c);
else
    PhaseID          = PhaseID_c;
end

if FigNo
    sfigure(FigNo); hold on
    tetramesh(EL2NOD(1:4,:)',GCOORD',...
              'EdgeColor','r',...
              'FaceColor','none',...
              'LineStyle','-','LineWidth',1);
end

end % END OF SUBFUNCTION tetmesh_refine_sph