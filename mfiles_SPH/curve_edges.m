function [GCOORD_curved,EL2NOD_curved,PointID_curved,GCOORD_SPH_curved,GCOORD_SPH_c,els_out_cone_no_cross_2pi,els_out_cone_cross_2pi,els_in_cone_no_iso,els_in_cone_iso] = ...
    curve_edges(GCOORD,EL2NOD,PointID,DB_indices,SETTINGS)
% Usage: [GCOORD_curved,EL2NOD_curved,PointID_curved,GCOORD_SPH_curved,GCOORD_SPH_c,els_out_cone_no_cross_2pi,els_out_cone_cross_2pi,els_in_cone_no_iso,els_in_cone_iso] = ...
%   curve_edges(GCOORD,EL2NOD,PointID,DB_indices,SETTINGS)
% 
% Purpose: 
%   Returns 10 nodel connectivity matrices (GCOORD in Cartesian coordinates
%   and GCOORD_SPH in spherical coordinates) after curving the element
%   edges.
%   
% Input:
%   GCOORD     : [matrix]    : Cartesian coordinates for a
%                              4-node/10-node/20-node tetrahedron mesh
%   EL2NOD     : [matrix]    : finite element connectivity matrix 
%                              (4 x nel)(10 x nel)(20 x nel)
%   PointID    : [vector]    : boundary indices for all nodes
%   DB_indices : [structure] : domain boundary indices
%   SETTINGS   : [structure] : model parameters
%
% Output:
%   GCOORD_curved             : [matrix] : Cartesian coordinates for a
%                                          10-node/20-node tetrahedron mesh
%                                          with curved edge elements
%   EL2NOD_curved             : [matrix] : finite element connectivity
%                                          matrix (10 x nel)(20 x nel)
%   PointID_curved            : [vector] : boundary indices for all nodes 
%   GCOORD_SPH_curved         : [matrix] : spherical coordinates for a
%                                          10-node/20-node tetrahedron mesh
%                                          with curved edge elements 
%   GCOORD_SPH_c              : [matrix] : spherical coordinates for a
%                                          4-node (vertices)
%   els_out_cone_no_cross_2pi : [vector] : elements outside the cone (4 
%                                          vertices outside the cone) and
%                                          not crossing phi = 2pi. These
%                                          elements have straight edges in
%                                          the original spherical frame
%   els_out_cone_cross_2pi    : [vector] : elements outside the cone (4 
%                                          vertices outside the cone) and
%                                          crossing phi = 2pi. These
%                                          elements have straight edges in
%                                          the rotated spherical frame
%                                          (180° around Z axis)
%   els_in_cone_no_iso        : [vector] : elements within the cone 
%                                          (4 vertices within the cone) +
%                                          elements crossing the cone
%                                          boundary and having only 1
%                                          vertex outside the cone 
%                                          (3 vertices within the cone +
%                                          1 vertex outside the cone) 
%                                          These elements have straight
%                                          edges in the rotated spherical
%                                          frame (90° around X axis). 
%   els_in_cone_iso           : [vector] : elements crossing the cone
%                                          boundary and having at least 1
%                                          edge (bar) outside the cone 
%                                          (2 vertices within the cone + 
%                                          2 vertices outside the cone OR 
%                                          1 vertex within the cone + 
%                                          3 vertices outside the cone).
%                                          These elements have curved edges
%                                          in the rotated spherical frame
%                                          (90° around X axis)
%
% JMT Jul 2016 : it needs to be coded for cubic curved edge elements

if nargin == 0
    pdir=pwd;cd('..');
    addpath([pwd '/SETUP_TEST/MESHES']);
    cd(pdir);
%     load('EarthShell_n562.mat')
    load('EarthShell_n4256.mat')
    GCOORD                = MESH.GCOORD;
    EL2NOD                = MESH.EL2NOD;
    PointID               = MESH.PointID;
    DB_indices{1}         = 301;
    DB_indices{2}         = 306;
    SETTINGS.element_type = 'quadratic';
    SETTINGS.theta_cone   = 45;
    SETTINGS.show_figs    = 1;
end

nnodel = size(EL2NOD,1);
if nnodel == 4
    return
elseif nnodel == 10 || 20 
    EL2NOD_c  = EL2NOD(1:4,:); % get the 4 nodel matrix
    nel       = max(max(EL2NOD_c));
    GCOORD_c  = GCOORD(:,1:nel); % 4 nodes per element mesh
    PointID_c = PointID(1:nel);
else
    error('incorrect nnodel')
end

% =========================================================================
% CONVERT 4 NODEL CARTESIAN MESH (x,y,z) WITH STRAIGHT EDGES TO 4 NODEL
% SPHERICAL MESH (theta,phi,r) WITH STRAIGHT EDGES
% =========================================================================
GCOORD_SPH_c = cartesian2spherical(GCOORD_c);

switch SETTINGS.element_type
    case 'quadratic'
        % =================================================================
        % CREATE A 10 NODEL MESH WITH CURVED EDGE ELEMENTS (CARTESIAN AND
        % SPHERICAL COORDINATES)
        % =================================================================
        [GCOORD_curved,GCOORD_SPH_curved,EL2NOD_curved,PointID_curved,els_out_cone_no_cross_2pi,els_out_cone_cross_2pi,els_in_cone_no_iso,els_in_cone_iso,~,~] = ...
            tetmesh_p1_to_p2_sph(GCOORD_c,GCOORD_SPH_c,EL2NOD_c,PointID_c,DB_indices,SETTINGS);
    case 'cubic'
        error('it needs to be coded')
        % =================================================================
        % CREATE A 20 NODEL SPHERICAL MESH (theta,phi,r) WITH CURVED EDGE
        % ELEMENTS
        % =================================================================
        [GCOORD_SPH_curved,EL2NOD_POL_curved,PointID_POL_curved] = ...
            trimesh_p1_to_p3_cyl(GCOORD_SPH_c,EL2NOD_c,PointID_c);
    otherwise
        error('SETTINGS.element_type must be "quadratic" or "cubic"')
end

if size(GCOORD_curved,2) < 600
    
    [els_cross_2pi,~] = check_phi(GCOORD_SPH_c,EL2NOD_c); % elements crossing phi = 2pi
    figure(61)
    clf
    subplot(1,2,1)
    plot_elements(GCOORD_c,EL2NOD_c,GCOORD,EL2NOD,els_cross_2pi)
    faces = [1 2 3 4];
    verts = [6371 0 6371; 0 0 6371; 0 0 -6371; 6371 0 -6371];
    patch('Faces',faces,'Vertices',verts,'FaceColor','w','FaceAlpha',0.1)
    title('Straight-edge elements crossing \phi = 2\pi')
    
    subplot(1,2,2)
    plot_elements(GCOORD_c,EL2NOD_c,GCOORD_curved,EL2NOD,els_cross_2pi)
    faces = [1 2 3 4];
    verts = [6371 0 6371; 0 0 6371; 0 0 -6371; 6371 0 -6371];
    patch('Faces',faces,'Vertices',verts,'FaceColor','w','FaceAlpha',0.1)
    title('Curved-edge elements crossing \phi = 2\pi')
    
end

if size(GCOORD,2) < 600
    figure(70)
    clf
    subplot(2,2,1)
    faceColor = [0.6875 0.8750 0.8984];
    tetramesh(EL2NOD_curved',GCOORD','FaceColor',faceColor,'FaceAlpha',0)
    hold on
    scatter3(GCOORD(1,:),GCOORD(2,:),GCOORD(3,:), ...
        'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
    axis equal
    view(142.5,30)
    grid on
    xlabel('x (km)')
    ylabel('y (km)')
    zlabel('z (km)')
    title('Straight edges Cartesian coordinates')
    
    subplot(2,2,3)
    [GCOORD_SPH_straight]    = cartesian2spherical(GCOORD);
    GCOORD_SPH_straight(3,:) = GCOORD_SPH_straight(3,:)/1000;
    tetramesh(EL2NOD_curved',GCOORD_SPH_straight','FaceColor',faceColor,'FaceAlpha',0)
    hold on
    scatter3(GCOORD_SPH_straight(1,:),GCOORD_SPH_straight(2,:),GCOORD_SPH_straight(3,:), ...
        'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
    axis([0 pi 0 2*pi 3 6.5])
    view(142.5,30)
    grid on
    xlabel('\theta (rad)')
    ylabel('\phi (rad)')
    zlabel('r (10^3 km)')
    title('Straight edges spherical coordinates')

    subplot(2,2,2)
    tetramesh(EL2NOD_curved',GCOORD_curved','FaceColor',faceColor,'FaceAlpha',0)
    hold on
    scatter3(GCOORD_curved(1,:),GCOORD_curved(2,:),GCOORD_curved(3,:), ...
        'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1])
    axis equal
    view(142.5,30)
    grid on
    xlabel('x (km)')
    ylabel('y (km)')
    zlabel('z (km)')
    title('Curved edges Cartesian coordinates')
    
    subplot(2,2,4)
    GCOORD_SPH_curved_plot      = GCOORD_SPH_curved;
    GCOORD_SPH_curved_plot(3,:) = GCOORD_SPH_curved_plot(3,:)/1000;
    tetramesh(EL2NOD_curved',GCOORD_SPH_curved_plot','FaceColor',faceColor,'FaceAlpha',0)
    hold on
    scatter3(GCOORD_SPH_curved_plot(1,:),GCOORD_SPH_curved_plot(2,:),GCOORD_SPH_curved_plot(3,:), ...
        'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1])
    axis([0 pi 0 2*pi 3 6.5])
    view(142.5,30)
    grid on
    xlabel('\theta (rad)')
    ylabel('\phi (rad)')
    zlabel('r (10^3 km)')
    title('Curved edges spherical coordinates')
end
end % END OF FUNCTION curved_edges

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function plot_elements(GCOORD_c,EL2NOD_c,GCOORD,EL2NOD,els)

GCOORD_els_vertices          = unique(GCOORD_c(:,EL2NOD_c(:,els))','rows','stable');
GCOORD_els_vertices          = GCOORD_els_vertices';
GCOORD_els                   = unique(GCOORD(:,EL2NOD(:,els))','rows','stable');
GCOORD_els                   = GCOORD_els';
lightGrey                    = 0.90*[1 1 1];
[x_sphere,y_sphere,z_sphere] = sphere(50);
x_sphere                     = x_sphere*6371;
y_sphere                     = y_sphere*6371;
z_sphere                     = z_sphere*6371;
faceColor                    = [0.6875 0.8750 0.8984];
axis equal
surface(x_sphere,y_sphere,z_sphere,'FaceColor','none','EdgeColor',lightGrey)
hold on
tetramesh(EL2NOD(:,els)',GCOORD','FaceColor',faceColor,'FaceAlpha',0.3)
axis([-6371 6371 -6371 6371 -6371 6371])
view(142.5,30)
grid on
scatter3(GCOORD_els(1,:),GCOORD_els(2,:),GCOORD_els(3,:), ...
    'MarkerEdgeColor','k','MarkerFaceColor',[0 1 0])
scatter3(GCOORD_els_vertices(1,:),GCOORD_els_vertices(2,:),GCOORD_els_vertices(3,:), ...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')

end % END OF SUBFUNCTION plot_elements