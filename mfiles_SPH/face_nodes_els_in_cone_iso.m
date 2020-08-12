function [GCOORD_face_nodes,EL2NOD14] = face_nodes_els_in_cone_iso(GCOORD,GCOORD_SPH_c,GCOORD_SPH,EL2NOD,els_in_cone_iso,SETTINGS)
% Usage: [GCOORD_face_nodes,EL2NOD14] = face_nodes_els_in_cone_iso(GCOORD,GCOORD_SPH_c,GCOORD_SPH,EL2NOD,els_in_cone_iso,SETTINGS)
%
% Purpose: 
%   Compute face nodes and connectivity for isoparametric elements. The
%   mid-face points are computed using the information of the 6 nodes (3
%   vertices + 3 midside points) that define each face.
%   Remember that an isoparametric element is defined here as an element
%   crossing the cone boundary and having at least 1 edge (bar) outside the
%   cone. There are 2 possible configurations for an isoparametric element: 
%       - 2 vertices within the cone + 2 vertices outside the cone
%       - 1 vertex within the cone + 3 vertices outside the cone
%   An isoparametric element has curved edges in the rotated spherical
%   frame (90° around X axis) since this element is merging both spherical
%   frames: rotated and original. 
%
% Input:
%   GCOORD            : [matrix]    : Cartesian coordinates of all nodes 
%                                     (4 nodel) (3 x nnod)
%   GCOORD_SPH_c      : [matrix]    : spherical coordinates of all nodes 
%                                     (4 nodel) (3 x nnod_c)
%   GCOORD_SPH        : [matrix]    : spherical coordinates of all nodes 
%                                     (10 nodel) (3 x nnod_c)
%   EL2NOD            : [matrix]    : connectivity matrix (10 x nel)
%   els_in_cone_iso   : [vector]    : isoparametric elements
%   SETTINGS          : [structure] : model parameters
%
% Output:
%   GCOORD_face_nodes : [matrix]    : Cartesian coordinates of all face
%                                     nodes  
%   EL2NOD14          : [matrix]    : connectivity matrix for isoparametric
%                                     elements (14 x nels_in_cone_iso)
%
% JMT Jul 2017
%

%=========================================================================================================================
% COMPUTE ALL FACE-NODES IN THE ORIGINAL FRAME
%=========================================================================================================================
EL2NOD_in_cone_iso = EL2NOD(:,els_in_cone_iso);  % EL2NOD for isoparametric elements
nnodel             = size(EL2NOD_in_cone_iso,1); % nodes per element
nel                = size(EL2NOD_in_cone_iso,2); % number of isoparametric elements
nnod0              = max(EL2NOD(:));             % total number of nodes
% Define faces using 6 nodes (according to sf_dsf_tet.m) since for
% isoparametric elements some edges are curved and some other are straight
faces = [2 3 4 6  7 10 ... % face 1
         1 3 4 9  7  8 ... % face 2
         1 2 4 5 10  8 ... % face 3
         1 2 3 5  6  9];   % face 4

FACE2NOD_all_faces = reshape( EL2NOD_in_cone_iso(faces',:),6,[] )';

% Find the faces that are shared by neighboring elements and return a unique list of faces.
[FACE2NOD_all_faces,~,ib] = unique_keep_order(FACE2NOD_all_faces);
nface                     = size(FACE2NOD_all_faces,1);

th_all_face_nodes = sum(reshape(GCOORD_SPH(1,FACE2NOD_all_faces'),6,nface))./6;
ph_all_face_nodes = sum(reshape(GCOORD_SPH(2,FACE2NOD_all_faces'),6,nface))./6;
r_all_face_nodes  = sum(reshape(GCOORD_SPH(3,FACE2NOD_all_faces'),6,nface))./6;

% element to bar connectivity after doubles have been merged (removed)
EL2FACE = reshape(uint32(ib),4,nel);

% write face node connectivity for isoparametric elements
EL2NOD14                      = EL2NOD_in_cone_iso;
EL2NOD14(nnodel+1:nnodel+4,:) = nnod0 + EL2FACE;

%======================================================================================================================================
% COMPUTE FACE-NODES THAT ARE OUTSIDE THE CONE BOUNDARY AND CROSSING PHI = 2PI IN THE 180° ROTATED FRAME AROUND Z AXIS 
% (WE USE ONLY THE 3 VERTICES TO COMPUTE THE FACE-NODES SINCE THESE FACES HAVE STRAIGHT EDGES IN THE 180° ROTATED FRAME AROUND Z AXIS) 
%======================================================================================================================================
TH2FACES          = reshape(GCOORD_SPH_c(1,FACE2NOD_all_faces(:,1:3)),[],3);
vertices_out_cone = sum(TH2FACES > SETTINGS.theta_cone*pi/180 & ...
                        TH2FACES < (180-SETTINGS.theta_cone)*pi/180,2); % vector for faces having x vertices outside the cone
                        % 1 --> one vertex is outside the cone boundary (this is not an isoparametric element)
                        % 2 --> two vertices are outside the cone boundary (no face outside the cone)
                        % 3 --> three vertices are outside the cone boundary (1 face outside the cone)
FACE2NOD_1_face_out_cone = FACE2NOD_all_faces(vertices_out_cone == 3,:);
PH2FACES          = reshape(GCOORD_SPH_c(2,FACE2NOD_1_face_out_cone(:,1:3)),[],3);
ang_dist          = zeros(3,size(FACE2NOD_1_face_out_cone,1)); 
ang_dist(1,:)     = abs(PH2FACES(:,1)-PH2FACES(:,2)); % angular distance between vertex 1 and vertex 2
ang_dist(2,:)     = abs(PH2FACES(:,2)-PH2FACES(:,3)); % angular distance between vertex 2 and vertex 3
ang_dist(3,:)     = abs(PH2FACES(:,3)-PH2FACES(:,1)); % angular distance between vertex 3 and vertex 1
ang_dist_longer_than_pi = ang_dist > pi; % angular distance longer than pi 
                                         % (means that some nodes are in the 1st/5th quadrant and some others nodes are in the 4th/8th quadrant)
faces_cross_2pi                    = sum(ang_dist_longer_than_pi,1) > 0; % boolean vector for elements crossing phi = 2pi
FACE2NOD_1_face_out_cone_cross_2pi = FACE2NOD_1_face_out_cone(faces_cross_2pi,:);

RR_Z_180_CCW                    = [-1  0  0 ; 0 -1  0 ; 0  0  1 ];     % 180° counterclockwise around the Z axis rotation matrix
GCOORD_rot_180                  = RR_Z_180_CCW * GCOORD;               % rotate coordinates
GCOORD_SPH_rot_180              = cartesian2spherical(GCOORD_rot_180); % convert to spherical coordinates
nface_out_cone_cross_2pi        = size(FACE2NOD_1_face_out_cone_cross_2pi,1);
th_faces_out_cone_cross_2pi_rot = sum(reshape(GCOORD_SPH_rot_180(1,FACE2NOD_1_face_out_cone_cross_2pi(:,1:3)'),3,nface_out_cone_cross_2pi))./3;
ph_faces_out_cone_cross_2pi_rot = sum(reshape(GCOORD_SPH_rot_180(2,FACE2NOD_1_face_out_cone_cross_2pi(:,1:3)'),3,nface_out_cone_cross_2pi))./3;
r_faces_out_cone_cross_2pi_rot  = sum(reshape(GCOORD_SPH_rot_180(3,FACE2NOD_1_face_out_cone_cross_2pi(:,1:3)'),3,nface_out_cone_cross_2pi))./3;

GCOORD_SPH_faces_out_cone_cross_2pi_rot = [th_faces_out_cone_cross_2pi_rot; ...
                                           ph_faces_out_cone_cross_2pi_rot; ...
                                           r_faces_out_cone_cross_2pi_rot];

GCOORD_faces_out_cone_cross_2pi_rot = spherical2cartesian(GCOORD_SPH_faces_out_cone_cross_2pi_rot); % convert to Cartesian coordinates
RR_Z_180_CW                         = [-1  0  0 ; 0 -1  0 ; 0  0  1 ];                     % 180° clockwise around the Z axis rotation matrix
GCOORD_faces_out_cone_cross_2pi     = RR_Z_180_CW * GCOORD_faces_out_cone_cross_2pi_rot;   % undo the rotation
GCOORD_SPH_faces_out_cone_cross_2pi = cartesian2spherical(GCOORD_faces_out_cone_cross_2pi); % convert to spherical coodinates

%=========================================================================================================================
% MERGE FACE NODES COMPUTED IN BOTH SPH FRAMES:
%   - FACE-NODES OUTSIDE THE CONE BOUNDARY (ORIGINAL FRAME)
%   - FACE-NODES OUTSIDE THE CONE BOUNDARY AND CROSSING PHI = 2PI (180° ROTATED FRAME)
%=========================================================================================================================
if ~isempty(FACE2NOD_1_face_out_cone_cross_2pi)
    % Compute positions of FACE2NOD_1_face_out_cone_cross_2pi in the FACE2NOD_all_faces list taking into account the possible permutations
    face_perms               = perms([3 2 1]);
    [~,LOCG]                 = ismember(FACE2NOD_1_face_out_cone_cross_2pi(:,face_perms(1,:)),FACE2NOD_all_faces(:,1:3),'rows');
    for i = 2:6
        [~,LOCG2]            = ismember(FACE2NOD_1_face_out_cone_cross_2pi(:,face_perms(i,:)),FACE2NOD_all_faces(:,1:3),'rows');
        LOCG(LOCG2 ~= 0)     = LOCG2(LOCG2 ~= 0);
    end
    th_all_face_nodes(LOCG) = GCOORD_SPH_faces_out_cone_cross_2pi(1,:);
    ph_all_face_nodes(LOCG) = GCOORD_SPH_faces_out_cone_cross_2pi(2,:);
    r_all_face_nodes(LOCG)  = GCOORD_SPH_faces_out_cone_cross_2pi(3,:);
end

%=========================================================================================================================
% COMPUTE FACE-NODES THAT ARE CROSSING THE CONE BOUNDARY IN THE 90° ROTATED FRAME AROUND X AXIS 
% (WE USE THE 3 VERTICES + 3 MIDPOINTS OF THE FACES CROSSING THE CONE BOUNDARY TO COMPUTE THEIR 
%  FACE-NODES SINCE THESE FACES HAVE CURVED EDGES IN THE 90° ROTATED FRAME AROUND X AXIS)  
%=========================================================================================================================
% Remove those ouside cone faces in the FACE2NOD_faces_crossing_cone list taking into account the possible permutations
face_perms           = perms([3 2 1]);
[~,LOCB]             = ismember(FACE2NOD_1_face_out_cone(:,face_perms(1,:)),FACE2NOD_all_faces(:,1:3),'rows');
for i = 2:6
    [~,LOCB2]        = ismember(FACE2NOD_1_face_out_cone(:,face_perms(i,:)),FACE2NOD_all_faces(:,1:3),'rows');
    LOCB(LOCB2 ~= 0) = LOCB2(LOCB2 ~= 0);
end
FACE2NOD_faces_crossing_cone = FACE2NOD_all_faces;
FACE2NOD_faces_crossing_cone(LOCB,:) = []; % remove those ouside cone faces in the FACE2NOD_faces_crossing_cone list

RR_X_90_CCW                = [ 1  0  0; 0  0 -1; 0  1  0 ];      % 90° counterclockwise around the X axis rotation matrix
GCOORD_rot_90              = RR_X_90_CCW * GCOORD;               % rotate coordinates
GCOORD_SPH_rot_90          = cartesian2spherical(GCOORD_rot_90); % convert to spherical coordinates
nface_crossing_cone        = size(FACE2NOD_faces_crossing_cone,1);
th_faces_crossing_cone_rot = sum(reshape(GCOORD_SPH_rot_90(1,FACE2NOD_faces_crossing_cone'),6,nface_crossing_cone))./6;
ph_faces_crossing_cone_rot = sum(reshape(GCOORD_SPH_rot_90(2,FACE2NOD_faces_crossing_cone'),6,nface_crossing_cone))./6;
r_faces_crossing_cone_rot  = sum(reshape(GCOORD_SPH_rot_90(3,FACE2NOD_faces_crossing_cone'),6,nface_crossing_cone))./6;

GCOORD_SPH_faces_crossing_cone_rot = [th_faces_crossing_cone_rot; ...
                                      ph_faces_crossing_cone_rot; ...
                                      r_faces_crossing_cone_rot];

GCOORD_faces_crossing_cone_rot = spherical2cartesian(GCOORD_SPH_faces_crossing_cone_rot); % convert to Cartesian coordinates
RR_X_90_CW                     = [ 1  0  0; 0  0  1; 0 -1  0 ];                   % 90° clockwise around the X axis rotation matrix
GCOORD_faces_crossing_cone     = RR_X_90_CW * GCOORD_faces_crossing_cone_rot;     % undo the rotation
GCOORD_SPH_faces_crossing_cone = cartesian2spherical(GCOORD_faces_crossing_cone); % convert to spherical coodinates

%=========================================================================================================================
% MERGE FACE NODES COMPUTED IN BOTH SPH FRAMES:
%   - FACE-NODES OUTSIDE THE CONE BOUNDARY (ORIGINAL FRAME)
%   - FACE-NODES CROSSING THE CONE BOUNDARY (90° ROTATED FRAME)
%=========================================================================================================================
if ~isempty(FACE2NOD_faces_crossing_cone)
    % Compute positions of FACE2NOD_faces_crossing_cone in the FACE2NOD_all_faces list taking into account the possible permutations
    face_perms               = perms([3 2 1]);
    [~,LOCF]                 = ismember(FACE2NOD_faces_crossing_cone(:,face_perms(1,:)),FACE2NOD_all_faces(:,1:3),'rows');
    for i = 2:6
        [~,LOCF2]            = ismember(FACE2NOD_faces_crossing_cone(:,face_perms(i,:)),FACE2NOD_all_faces(:,1:3),'rows');
        LOCF(LOCF2 ~= 0)     = LOCF2(LOCF2 ~= 0);
    end
    th_all_face_nodes(LOCF) = GCOORD_SPH_faces_crossing_cone(1,:);
    ph_all_face_nodes(LOCF) = GCOORD_SPH_faces_crossing_cone(2,:);
    r_all_face_nodes(LOCF)  = GCOORD_SPH_faces_crossing_cone(3,:);
end

GCOORD_SPH_face_nodes = [th_all_face_nodes; ph_all_face_nodes; r_all_face_nodes];
GCOORD_face_nodes     = spherical2cartesian(GCOORD_SPH_face_nodes); % convert to Cartesian coordinates

% %=========================================================================================================================
% % PLOT
% %=========================================================================================================================
% figure(742);clf
% GCOORD_plot = [GCOORD GCOORD_face_nodes];
% els_plot    = 115;%2; %13; %75;
% lightMag    = 0.80*[1 0.7 1]; % colour for the cone boundary
% lightGrey   = 0.90*[1 1 1];   % colour for the shell boundaries
% lightWhi    = [1 1 1];       % colour for elements
% th          = SETTINGS.theta_cone*pi/180;
% h           = 6371;
% r           = h*tan(th);
% [R,A]       = meshgrid(linspace(0,r,2),linspace(0,2*pi,21));
% X           = R .* cos(A);
% Y           = R .* sin(A);
% Z           = R/tan(th);
% surface(X,Y,Z,'FaceColor','none','EdgeColor',lightMag)
% hold on
% surface(X,Y,-Z,'FaceColor','none','EdgeColor',lightMag)
% % Plot the shell
% [x_sph,y_sph,z_sph] = sphere(20);
% x_sph               = x_sph*3471;
% y_sph               = y_sph*3471;
% z_sph               = z_sph*3471;
% axis equal
% surface(x_sph,y_sph,z_sph,'FaceColor','none','EdgeColor',lightGrey)
% [x_sph,y_sph,z_sph] = sphere(30);
% x_sph               = x_sph*6371;
% y_sph               = y_sph*6371;
% z_sph               = z_sph*6371;
% surface(x_sph,y_sph,z_sph,'FaceColor','none','EdgeColor',lightGrey)
% axis([-6371 6371 -6371 6371 -6371 6371])
% view(142.5,30)
% grid on
% hold on
% tetramesh(EL2NOD14(1:4,els_plot)',GCOORD','FaceColor',lightWhi,'FaceAlpha',0.3)
% color_local_vertex_1 = [1 0 0]; % red
% color_local_vertex_2 = [1 1 0]; % yellow
% color_local_vertex_3 = [0 1 0]; % green
% color_local_vertex_4 = [0 0 1]; % blue
% scatter3(GCOORD(1,EL2NOD14(1,els_plot)), ...
%          GCOORD(2,EL2NOD14(1,els_plot)), ...
%          GCOORD(3,EL2NOD14(1,els_plot)),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_1)
% scatter3(GCOORD(1,EL2NOD14(2,els_plot)), ...
%          GCOORD(2,EL2NOD14(2,els_plot)), ...
%          GCOORD(3,EL2NOD14(2,els_plot)),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_2)
% scatter3(GCOORD(1,EL2NOD14(3,els_plot)), ...
%          GCOORD(2,EL2NOD14(3,els_plot)), ...
%          GCOORD(3,EL2NOD14(3,els_plot)),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_3)
% scatter3(GCOORD(1,EL2NOD14(4,els_plot)), ...
%          GCOORD(2,EL2NOD14(4,els_plot)), ...
%          GCOORD(3,EL2NOD14(4,els_plot)),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_4)
% 
% scatter3(GCOORD_plot(1,EL2NOD14(11,els_plot)), ...
%          GCOORD_plot(2,EL2NOD14(11,els_plot)), ...
%          GCOORD_plot(3,EL2NOD14(11,els_plot)),'MarkerEdgeColor',color_local_vertex_1,'MarkerFaceColor',color_local_vertex_1)
% scatter3(GCOORD_plot(1,EL2NOD14(12,els_plot)), ...
%          GCOORD_plot(2,EL2NOD14(12,els_plot)), ...
%          GCOORD_plot(3,EL2NOD14(12,els_plot)),'MarkerEdgeColor',color_local_vertex_2,'MarkerFaceColor',color_local_vertex_2)
% scatter3(GCOORD_plot(1,EL2NOD14(13,els_plot)), ...
%          GCOORD_plot(2,EL2NOD14(13,els_plot)), ...
%          GCOORD_plot(3,EL2NOD14(13,els_plot)),'MarkerEdgeColor',color_local_vertex_3,'MarkerFaceColor',color_local_vertex_3)
% scatter3(GCOORD_plot(1,EL2NOD14(14,els_plot)), ...
%          GCOORD_plot(2,EL2NOD14(14,els_plot)), ...
%          GCOORD_plot(3,EL2NOD14(14,els_plot)),'MarkerEdgeColor',color_local_vertex_4,'MarkerFaceColor',color_local_vertex_4)

end % END OF FUNCTION face_nodes_els_in_cone_iso