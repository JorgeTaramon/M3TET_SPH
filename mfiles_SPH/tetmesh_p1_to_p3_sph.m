function [GCOORD,GCOORD_SPH,EL2NOD] = ...
    tetmesh_p1_to_p3_sph(GCOORD4,GCOORD_SPH4,EL2NOD4,r_cmb,r_surf,SETTINGS)
% Usage: [GCOORD,GCOORD_SPH,EL2NOD] = ...
%   tetmesh_p1_to_p3_sph(GCOORD4,GCOORD_SPH4,EL2NOD4,r_cmb,r_surf,SETTINGS)
%
% Purpose:
%   Creates a cubic (20-node) tetrahedra mesh from a linear (4-node)
%   tetrahedra mesh in spherical coordinates
%
% Input:
%   GCOORD4     : [matrix]    : Cartesian coordinates of all nodes 
%                               (4 nodel) (3 x nnod_c)
%   GCOORD_SPH4 : [matrix]    : spherical coordinates of all nodes 
%                               (4 nodel) (3 x nnod_c)
%   EL2NOD4     : [matrix]    : connectivity matrix (4 x nel)
%   r_cmb       : [scalar]    : radius for CMB (km)
%   r_surf      : [scalar]    : radius for surface (km)
%   SETTINGS    : [structure] : model parameters
%
% Output:
%   GCOORD             : [matrix]   : Cartesian coordinates of all
%                                            nodes new mesh (3 x ?)
%   GCOORD_SPH                : [matrix]   : spherical coordinates of all
%                                            nodes new mesh (3 x ?)
%   EL2NOD                    : [matrix]   : connectivity matrix (10 x nel)
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH March 2013
% JMT Dec 2016: Compute curved edges for a spherical mesh
%

%======================================================================================================
% COMPUTE ELEMENTS IN RELATION WITH THE CONE AND CROSSING PHI = 2PI
%======================================================================================================
[els_in_cone,els_within_cone,els_cone_bnd,els_out_cone,els_in_cone_iso] = check_els_in_cone(GCOORD_SPH4,EL2NOD4,SETTINGS.theta_cone);
[els_cross_2pi,~]         = check_phi(GCOORD_SPH4,EL2NOD4(:,els_out_cone));
els_out_cone_cross_2pi    = els_out_cone(els_cross_2pi);
els_out_cone_no_cross_2pi = els_out_cone;
els_out_cone_no_cross_2pi(ismember(els_out_cone_no_cross_2pi,els_out_cone_cross_2pi)) = [];
els_in_cone_no_iso        = els_in_cone;
els_in_cone_no_iso(ismember(els_in_cone_no_iso,els_in_cone_iso)) = [];

ELS.els_out_cone_no_cross_2pi = els_out_cone_no_cross_2pi;
ELS.els_out_cone_cross_2pi    = els_out_cone_cross_2pi;
ELS.els_in_cone_no_iso        = els_in_cone_no_iso;
ELS.els_in_cone_iso           = els_in_cone_iso;

%======================================================================================================
% COMPUTE EDGE NODES AND FACE NODES IN SPHERICAL COORDINATES FOR ELEMENTS OUTSIDE 
% THE CONE AND CROSSING PHI = 2PI AFTER A 180° ROTATION AROUND Z AXIS 
%======================================================================================================
if ~isempty(els_out_cone_cross_2pi)
    [GCOORD_SPH_bars_els_cross_2pi,EDGE2NOD_els_cross_2pi,GCOORD_SPH_face_nodes_els_cross_2pi,FACE2NOD_els_cross_2pi] = ...
        compute_cubic_nodes_sph_rot_frame_els_cross_2pi(GCOORD4,EL2NOD4,els_out_cone_cross_2pi);
else
    GCOORD_SPH_bars_els_cross_2pi       = [];
    EDGE2NOD_els_cross_2pi              = [];
    GCOORD_SPH_face_nodes_els_cross_2pi = [];
    FACE2NOD_els_cross_2pi              = [];
end

%======================================================================================================
% COMPUTE EDGE NODES AND FACE NODES IN SPHERICAL COORDINATES FOR ELEMENTS IN THE CONE 
% (WITHIN THE CONE BOUNDARY + CROSSING THE CONE BOUNDARY) AFTER A 90° ROTATION AROUND X AXIS 
%======================================================================================================
if ~isempty(els_in_cone)
    [GCOORD_SPH_bars_els_in_cone,EDGE2NOD_els_in_cone,GCOORD_SPH_face_nodes_els_in_cone,FACE2NOD_els_in_cone] = ...
        compute_cubic_nodes_sph_rot_frame_els_in_cone(GCOORD4,GCOORD_SPH4,EL2NOD4,els_in_cone,els_in_cone_iso,SETTINGS.theta_cone);
else
    GCOORD_SPH_bars_els_in_cone       = [];
    EDGE2NOD_els_in_cone              = [];
    GCOORD_SPH_face_nodes_els_in_cone = [];
    FACE2NOD_els_in_cone              = [];
end

%======================================================================================================
% COMPUTE EDGE NODES AND FACE NODES IN SPHERICAL COORDINATES (CURVED EDGES) FOR ALL ELEMENTS
%======================================================================================================
[GCOORD,GCOORD_SPH,EL2NOD] = compute_cubic_nodes_sph_rot_frame_all_els...
    (GCOORD_SPH4,EL2NOD4,ELS,GCOORD_SPH_bars_els_in_cone,EDGE2NOD_els_in_cone,GCOORD_SPH_face_nodes_els_in_cone,FACE2NOD_els_in_cone,...
     GCOORD_SPH_bars_els_cross_2pi,EDGE2NOD_els_cross_2pi,GCOORD_SPH_face_nodes_els_cross_2pi,FACE2NOD_els_cross_2pi);

% Fig_num    = 1;
% els        = [1 2 3];
% lightCyan  = 0.8*[0 1 1]; % colour for the elements within the cone and isoparametric
% lightRed   = [0.8 0 0];   % colour for nodes
% theta_cone = 45;
% plot_points_and_elements(GCOORD,GCOORD_SPH,EL2NOD,GCOORD,GCOORD_SPH,els,theta_cone,lightCyan,lightRed,Fig_num)

% % %======================================================================================================
% % % DATA FOR OUTPUT
% % %======================================================================================================
% % % Calculate boundary index for the new nodes
% % FigNo   = 0;
% % PointID = calc_PointIDs(GCOORD,EL2NOD,PointID_c,DB_indices,FigNo);
% % 
% % % Make sure that boundary nodes are exactly on the boundaries (the radius of boundary nodes can shift
% % % slightly (1e-8) from the actual boundary because of transforming between Cartesian and spherical coord)
% % GCOORD_SPH(3,PointID==301) = r_cmb;
% % GCOORD_SPH(3,PointID==306) = r_surf;

end % END OF FUNCTION tetmesh_p1_to_p3_sph

% #####################################################################################################
%                                           SUB-FUNCTIONS
% #####################################################################################################

function phiBarEnds = check_ang_dist_phi(GCOORD_SPH_c,bar2node,phiBarEnds)
% Usage: phiBarEnds = check_ang_dist_phi(GCOORD_SPH_c,bar2node,phiBarEnds)
%
% Purpose:
%   Find which bars have phi angles (of the end-nodes) separated an angular 
%   distance bigger than pi. 
%   Then, change those values of phi < pi by phi + 2pi
%
% Input:
%   GCOORD_SPH_c : [matrix] : spherical coordinates of the 4-node
%                             tetrahedron mesh
%   bar2node     : [matrix] : pointer bars defined by their end-nodes
%   phiBarEnds   : [matrix] : phi angles for the end-nodes of bars
%
% Output:
%   phiBarEnds   : [matrix] : phi angles for the end-nodes of bars
%
% JMT Jul 2016

ang_dist = abs(GCOORD_SPH_c(2,bar2node(:,1))-GCOORD_SPH_c(2,bar2node(:,2)));
ind      = find(ang_dist > pi); % indices for those bars having the phi angles (of the ends) separated more than pi
for i=1:size(ind,2)
    % Change those values of phi < pi by phi + 2pi
    phiBarEnds(phiBarEnds(:,ind(i))<pi,ind(i)) = phiBarEnds(phiBarEnds(:,ind(i))<pi,ind(i)) + 2*pi;
end
end % END OF SUBFUNCTION check_ang_dist_phi

% #####################################################################################################

function [GCOORD_SPH_bars_els_cross_2pi,EDGE2NOD_els_cross_2pi,GCOORD_SPH_face_nodes_els_cross_2pi,FACE2NOD_els_cross_2pi] = ...
    compute_cubic_nodes_sph_rot_frame_els_cross_2pi(GCOORD4,EL2NOD4,els_out_cone_cross_2pi)
% Usage: [GCOORD_SPH_bars_els_cross_2pi,EDGE2NOD_els_cross_2pi,GCOORD_SPH_face_nodes_els_cross_2pi,FACE2NOD_els_cross_2pi] = ...
%   compute_cubic_nodes_sph_rot_frame_els_cross_2pi(GCOORD4,EL2NOD4,els_out_cone_cross_2pi)
%
% Purpose: 
%   
%
% Input:
%   
%
% Output:
%  
%
% JMT Dec 2016

%=========================================================================================================================
% COMPUTE NODES AT 1/3 and 2/3 IN THE EDGES
%=========================================================================================================================
% COMPUTE EDGES FOR ALL ELEMENTS CROSSING PHI = 2PI
% ------------------------------------------------------------------------------------------------------------------------
EL2NOD4_cross_2pi            = EL2NOD4(:,els_out_cone_cross_2pi); % connectivity matrix for elements crossing phi = 2pi
% Create a pointer edges defined by their end-nodes (e.g. bar 1 has end-nodes [1 2], bar 2 has [2 3],...)
edges                        = [ 1  2 ... % 1st edge of parent
                                 2  3 ... % 2nd edge of parent
                                 3  4 ... % 3rd edge of parent
                                 4  1 ... % 4th edge of parent
                                 1  3 ... % 5th edge of parent
                                 4  2 ];  % 6th edge of parent
EDGE2NOD_els_cross_2pi       = reshape( EL2NOD4_cross_2pi(edges',:),2,[] )';
[EDGE2NOD_els_cross_2pi,~,~] = unique_keep_order(EDGE2NOD_els_cross_2pi); % find the edges that are shared by neighboring elements and return a unique list of edges
nedge                        = size(EDGE2NOD_els_cross_2pi,1); % number of unique edges (i.e. number of new edge nodes that have to be generated)

% COMPUTE COORDINATES FOR NEW NODES ON THE EDGES
% ----------------------------------------------
RR_Z_180_CCW    = [-1  0  0;   0 -1  0; 0  0  1 ];  % 180° counterclockwise around the Z axis rotation matrix
GCOORD4_rot     = RR_Z_180_CCW * GCOORD4;           % rotate coordinates
GCOORD_SPH4_rot = cartesian2spherical(GCOORD4_rot); % convert to spherical coordinates

% Coordinates of the nodes defining each bar's end points
thBarEnds_els_cross_2pi   = reshape(GCOORD_SPH4_rot(1,EDGE2NOD_els_cross_2pi'),2,nedge);
phBarEnds_els_cross_2pi   = reshape(GCOORD_SPH4_rot(2,EDGE2NOD_els_cross_2pi'),2,nedge);
rBarEnds_els_cross_2pi    = reshape(GCOORD_SPH4_rot(3,EDGE2NOD_els_cross_2pi'),2,nedge);

% Create two new nodes on each edge
LthBar_els_cross_2pi      = thBarEnds_els_cross_2pi(2,:) - thBarEnds_els_cross_2pi(1,:);
thBarThirds_els_cross_2pi = [thBarEnds_els_cross_2pi(1,:)+(1/3)*LthBar_els_cross_2pi
                             thBarEnds_els_cross_2pi(1,:)+(2/3)*LthBar_els_cross_2pi];
thBarThirds_els_cross_2pi = thBarThirds_els_cross_2pi(:)';
LphBar_els_cross_2pi      = phBarEnds_els_cross_2pi(2,:) - phBarEnds_els_cross_2pi(1,:);
phBarThirds_els_cross_2pi = [phBarEnds_els_cross_2pi(1,:)+(1/3)*LphBar_els_cross_2pi
                             phBarEnds_els_cross_2pi(1,:)+(2/3)*LphBar_els_cross_2pi];
phBarThirds_els_cross_2pi = phBarThirds_els_cross_2pi(:)';
LrBar_els_cross_2pi       = rBarEnds_els_cross_2pi(2,:) - rBarEnds_els_cross_2pi(1,:);
rBarThirds_els_cross_2pi  = [rBarEnds_els_cross_2pi(1,:)+(1/3)*LrBar_els_cross_2pi
                             rBarEnds_els_cross_2pi(1,:)+(2/3)*LrBar_els_cross_2pi];
rBarThirds_els_cross_2pi  = rBarThirds_els_cross_2pi(:)';

GCOORD_SPH_bars_els_cross_2pi_rot = [thBarThirds_els_cross_2pi; ...
                                     phBarThirds_els_cross_2pi; ...
                                     rBarThirds_els_cross_2pi];

GCOORD_bars_els_cross_2pi_rot = spherical2cartesian(GCOORD_SPH_bars_els_cross_2pi_rot); % convert to Cartesian coordinates
RR_Z_180_CW                   = [-1  0  0;   0 -1  0; 0  0  1 ];  % 180° clockwise around the Z axis rotation matrix
GCOORD_bars_els_cross_2pi     = RR_Z_180_CW * GCOORD_bars_els_cross_2pi_rot;     % undo the rotation
GCOORD_SPH_bars_els_cross_2pi = cartesian2spherical(GCOORD_bars_els_cross_2pi); % convert to spherical coodinates

%=========================================================================================================================
% COMPUTE FACE NODES
%=========================================================================================================================
% COMPUTE FACES FOR ALL ELEMENTS CROSSING PHI = 2PI
% ------------------------------------------------------------------------------------------------------------------------
[GCOORD_SPH_face_nodes_els_cross_2pi_rot,~,FACE2NOD_els_cross_2pi]   = tetmesh_add_facenods(GCOORD_SPH4_rot,EL2NOD4_cross_2pi);
GCOORD_SPH_face_nodes_els_cross_2pi_rot(:,1:size(GCOORD_SPH4_rot,2)) = []; % select only coordinates for face nodes

GCOORD_face_nodes_els_cross_2pi_rot = spherical2cartesian(GCOORD_SPH_face_nodes_els_cross_2pi_rot); % convert to Cartesian coordinates
GCOORD_face_nodes_els_cross_2pi     = RR_Z_180_CW * GCOORD_face_nodes_els_cross_2pi_rot;            % undo the rotation
GCOORD_SPH_face_nodes_els_cross_2pi = cartesian2spherical(GCOORD_face_nodes_els_cross_2pi);         % convert to spherical coodinates

end % END OF SUBFUNCTION mid_side_nodes_sph_rot_frame_els_cross_2pi

% #####################################################################################################

function [GCOORD_SPH_bars_els_in_cone,EDGE2NOD_els_in_cone,GCOORD_SPH_face_nodes_els_in_cone,FACE2NOD_els_in_cone] = ...
    compute_cubic_nodes_sph_rot_frame_els_in_cone(GCOORD4,GCOORD_SPH4,EL2NOD4,els_in_cone,els_in_cone_iso,theta_cone)
% Usage: [GCOORD_SPH_bars_els_in_cone,EDGE2NOD_in_cone,GCOORD_SPH_face_nodes_els_in_cone,FACE2NOD_in_cone] = ...
%   compute_cubic_nodes_sph_rot_frame_els_in_cone(GCOORD4,GCOORD_SPH4,EL2NOD4,els_in_cone,els_in_cone_iso,theta_cone)
%
% Purpose: 
%   
%
% Input:
%   
%
% Output:
%  
%
% JMT Dec 2016

%=========================================================================================================================
% COMPUTE NODES AT 1/3 and 2/3 IN THE EDGES
%=========================================================================================================================
% COMPUTE EDGES FOR ALL ELEMENTS INSIDE THE CONE 
% ------------------------------------------------------------------------------------------------------------------------
EL2NOD4_in_cone            = EL2NOD4(:,els_in_cone); % connectivity matrix for elements inside the cone
% Create a pointer edges defined by their end-nodes (e.g. bar 1 has end-nodes [1 2], bar 2 has [2 3],...)
edges                      = [ 1  2 ... % 1st edge of parent
                               2  3 ... % 2nd edge of parent
                               3  4 ... % 3rd edge of parent
                               4  1 ... % 4th edge of parent
                               1  3 ... % 5th edge of parent
                               4  2 ];  % 6th edge of parent
EDGE2NOD_els_in_cone       = reshape( EL2NOD4_in_cone(edges',:),2,[] )';
[EDGE2NOD_els_in_cone,~,~] = unique_keep_order(EDGE2NOD_els_in_cone); % find the edges that are shared by neighboring elements and return a unique list of edges

if ~isempty(els_in_cone_iso)
    % COMPUTE EDGES FOR ALL ISOPARAMETRIC ELEMENTS (ELEMENTS INSIDE THE CONE HAVING AT LEAST 1 EDGE OUTSIDE THE CONE BOUNDARY)
    % ------------------------------------------------------------------------------------------------------------------------
    EL2NOD4_in_cone_iso        = EL2NOD4(:,els_in_cone_iso);
    EDGE2NOD_in_cone_iso       = reshape( EL2NOD4_in_cone_iso(edges',:),2,[] )';
    [EDGE2NOD_in_cone_iso,~,~] = unique_keep_order(EDGE2NOD_in_cone_iso);
    
    % Select those edges outside the cone boundary (both edge ends are outside of the cone boundary)
    TH2EDGES                   = reshape(GCOORD_SPH4(1,EDGE2NOD_in_cone_iso),[],2);
    theta_bars_out_cone        = sum(TH2EDGES > theta_cone*pi/180 & TH2EDGES < (180-theta_cone)*pi/180,2); % vector for edges having both ends outside the cone boundary
                                    % 0 --> both ends bar are inside the cone boundary
                                    % 1 --> one end bar is inside the cone boundary and the other one is outside the cone boundary
                                    % 2 --> both ends bar are outside the cone boundary
    EDGE2NOD_in_cone_iso_with_2_bars_out_cone = EDGE2NOD_in_cone_iso(theta_bars_out_cone == 2,:);
    
    % Remove those ouside cone edges in the EDGE2NOD_in_cone list
    [~,LOCA]        = ismember(EDGE2NOD_in_cone_iso_with_2_bars_out_cone,EDGE2NOD_els_in_cone,'rows');
    [~,LOCA2]       = ismember(fliplr(EDGE2NOD_in_cone_iso_with_2_bars_out_cone),EDGE2NOD_els_in_cone,'rows'); % in case some edges were flipped when using unique_keep_order
    LOCA(LOCA == 0) = LOCA2(LOCA2 ~= 0); clear LOCA2;
    EDGE2NOD_els_in_cone(LOCA,:) = []; % remove those ouside cone edges in the EDGE2NOD_in_cone list
end

% COMPUTE COORDINATES FOR NEW NODES ON THE EDGES
% ----------------------------------------------
RR_X_90_CCW     = [ 1  0  0; 0  0 -1; 0  1  0 ];    % 90° counterclockwise around the X axis rotation matrix
GCOORD4_rot     = RR_X_90_CCW * GCOORD4;            % rotate coordinates
GCOORD_SPH4_rot = cartesian2spherical(GCOORD4_rot); % convert to spherical coordinates

% Coordinates of the nodes defining each bar's end points
nedge                   = size(EDGE2NOD_els_in_cone,1); % number of unique edges (i.e. number of new edge nodes that have to be generated)
thBarEnds_els_in_cone   = reshape(GCOORD_SPH4_rot(1,EDGE2NOD_els_in_cone'),2,nedge);
phBarEnds_els_in_cone   = reshape(GCOORD_SPH4_rot(2,EDGE2NOD_els_in_cone'),2,nedge);
rBarEnds_els_in_cone    = reshape(GCOORD_SPH4_rot(3,EDGE2NOD_els_in_cone'),2,nedge);

% Create two new nodes on each edge
LthBar_els_in_cone      = thBarEnds_els_in_cone(2,:) - thBarEnds_els_in_cone(1,:);
thBarThirds_els_in_cone = [thBarEnds_els_in_cone(1,:)+(1/3)*LthBar_els_in_cone
                           thBarEnds_els_in_cone(1,:)+(2/3)*LthBar_els_in_cone];
thBarThirds_els_in_cone = thBarThirds_els_in_cone(:)';
LphBar_els_in_cone      = phBarEnds_els_in_cone(2,:) - phBarEnds_els_in_cone(1,:);
phBarThirds_els_in_cone = [phBarEnds_els_in_cone(1,:)+(1/3)*LphBar_els_in_cone
                           phBarEnds_els_in_cone(1,:)+(2/3)*LphBar_els_in_cone];
phBarThirds_els_in_cone = phBarThirds_els_in_cone(:)';
LrBar_els_in_cone       = rBarEnds_els_in_cone(2,:) - rBarEnds_els_in_cone(1,:);
rBarThirds_els_in_cone  = [rBarEnds_els_in_cone(1,:)+(1/3)*LrBar_els_in_cone
                           rBarEnds_els_in_cone(1,:)+(2/3)*LrBar_els_in_cone];
rBarThirds_els_in_cone  = rBarThirds_els_in_cone(:)';

GCOORD_SPH_bars_els_in_cone_rot = [thBarThirds_els_in_cone; ...
                                   phBarThirds_els_in_cone; ...
                                   rBarThirds_els_in_cone];

GCOORD_bars_els_in_cone_rot = spherical2cartesian(GCOORD_SPH_bars_els_in_cone_rot); % convert to Cartesian coordinates
RR_X_90_CW                  = [ 1  0  0; 0  0  1; 0 -1  0 ]; % 90° clockwise around the X axis rotation matrix
GCOORD_bars_els_in_cone     = RR_X_90_CW * GCOORD_bars_els_in_cone_rot; % undo the rotation
GCOORD_SPH_bars_els_in_cone = cartesian2spherical(GCOORD_bars_els_in_cone); % convert to spherical coodinates


%=========================================================================================================================
% COMPUTE FACE NODES
%=========================================================================================================================
% COMPUTE FACES FOR ALL ELEMENTS INSIDE THE CONE 
% ------------------------------------------------------------------------------------------------------------------------
[~,~,FACE2NOD_els_in_cone] = tetmesh_add_facenods(GCOORD_SPH4_rot,EL2NOD4_in_cone);
% GCOORD_SPH_face_els_in_cone_rot(:,1:size(GCOORD_SPH4_rot,2)) = [];  % remove coordinates for vertices (keep just the coordinates for the face nodes)
% EL2FACE_els_in_cone = repmat([1 2 3 4]',1,size(EL2NOD4_in_cone,2)); % face id's for each element inside the cone

if ~isempty(els_in_cone_iso)
    % COMPUTE FACES FOR ALL ISOPARAMETRIC ELEMENTS (ELEMENTS INSIDE THE CONE HAVING AT LEAST 1 EDGE OUTSIDE THE CONE BOUNDARY)
    % ------------------------------------------------------------------------------------------------------------------------
    % Select those faces outside the cone boundary (three vertices defining each face are outside of the cone boundary)
    [~,~,FACE2NOD_els_in_cone_iso] = tetmesh_add_facenods(GCOORD_SPH4_rot,EL2NOD4_in_cone_iso);
    TH2FACES                       = reshape(GCOORD_SPH4(1,FACE2NOD_els_in_cone_iso),[],3);
    theta_faces_out_cone           = sum(TH2FACES > theta_cone*pi/180 & TH2FACES < (180-theta_cone)*pi/180,2); % vector for faces having the 3 vertices outside the cone
                                        % 1 --> one vertex is outside the cone boundary
                                        % 2 --> two vertices are outside the cone boundary
                                        % 3 --> three vertices are outside the cone boundary
    FACE2NOD_els_in_cone_iso_with_1_face_out_cone = FACE2NOD_els_in_cone_iso(theta_faces_out_cone == 3,:);
    
    % Remove those ouside cone faces in the FACE2NOD_in_cone list taking into account the possible permutations
    face_perms               = perms([3 2 1]);
    [~,LOCB]                 = ismember(FACE2NOD_els_in_cone_iso_with_1_face_out_cone(:,face_perms(1,:)),FACE2NOD_els_in_cone,'rows');
    for i = 2:6
        [~,LOCB2]            = ismember(FACE2NOD_els_in_cone_iso_with_1_face_out_cone(:,face_perms(i,:)),FACE2NOD_els_in_cone,'rows');
        LOCB(LOCB2 ~= 0)     = LOCB2(LOCB2 ~= 0);
    end
    FACE2NOD_els_in_cone(LOCB,:) = []; % remove those ouside cone faces in the FACE2NOD_in_cone list
end

% COMPUTE COORDINATES FOR FACE NODES
% ----------------------------------
nface                     = size(FACE2NOD_els_in_cone,1);
th_face_nodes_els_in_cone = sum(reshape(GCOORD_SPH4_rot(1,FACE2NOD_els_in_cone'),3,nface))./3;
ph_face_nodes_els_in_cone = sum(reshape(GCOORD_SPH4_rot(2,FACE2NOD_els_in_cone'),3,nface))./3;
r_face_nodes_els_in_cone  = sum(reshape(GCOORD_SPH4_rot(3,FACE2NOD_els_in_cone'),3,nface))./3;

GCOORD_SPH_face_nodes_els_in_cone_rot = [th_face_nodes_els_in_cone; ...
                                         ph_face_nodes_els_in_cone; ...
                                         r_face_nodes_els_in_cone];

GCOORD_face_nodes_els_in_cone_rot = spherical2cartesian(GCOORD_SPH_face_nodes_els_in_cone_rot); % convert to Cartesian coordinates
GCOORD_face_nodes_els_in_cone     = RR_X_90_CW * GCOORD_face_nodes_els_in_cone_rot; % undo the rotation
GCOORD_SPH_face_nodes_els_in_cone = cartesian2spherical(GCOORD_face_nodes_els_in_cone); % convert to spherical coodinates

end % END OF SUBFUNCTION compute_cubic_nodes_sph_rot_frame_els_in_cone

% #####################################################################################################

function [GCOORD,GCOORD_SPH,EL2NOD] = compute_cubic_nodes_sph_rot_frame_all_els...
    (GCOORD_SPH4,EL2NOD4,ELS,GCOORD_SPH_bars_els_in_cone,EDGE2NOD_in_cone,GCOORD_SPH_face_nodes_els_in_cone,FACE2NOD_in_cone,...
     GCOORD_SPH_bars_els_cross_2pi,EDGE2NOD_els_cross_2pi,GCOORD_SPH_face_nodes_els_cross_2pi,FACE2NOD_els_cross_2pi)
% Usage: [GCOORD_curved,GCOORD_SPH,EL2NOD] = mid_side_nodes_sph_frame_all_els...
%   (GCOORD_SPH4,EL2NOD4,ELS,GCOORD_SPH_bars_els_in_cone,EDGE2NOD_in_cone,GCOORD_SPH_face_nodes_els_in_cone,FACE2NOD_in_cone,...
%    GCOORD_SPH_bars_els_cross_2pi,EDGE2NOD_els_cross_2pi,GCOORD_SPH_face_nodes_els_cross_2pi,FACE2NOD_els_cross_2pi)
%
% Purpose: 
%  
%
% Input:
%  
% Output:
%
% JMT Dec 2016

%=========================================================================================================================
% COMPUTE NODES AT 1/3 and 2/3 IN THE EDGES
%=========================================================================================================================
% Size of coarse mesh
nel              = size(EL2NOD4,2);     % number of elements
nnod4            = size(GCOORD_SPH4,2); % number of nodes
% Create a pointer edges defined by their end-nodes (e.g. bar 1 has end-nodes [1 2], bar 2 has [2 3],...)
edges                  = [ 1  2 ... % 1st edge of parent
                           2  3 ... % 2nd edge of parent
                           3  4 ... % 3rd edge of parent
                           4  1 ... % 4th edge of parent
                           1  3 ... % 5th edge of parent
                           4  2 ];  % 6th edge of parent
EDGE2NOD        = reshape( EL2NOD4(edges',:),2,[] )';
% Find the edges that are shared by neighboring elements and return a unique list of edges.
[EDGE2NOD,~,ib] = unique_keep_order(EDGE2NOD);
nedge           = size(EDGE2NOD,1); % number of unique edges (i.e. number of new edge nodes that have to be generated)

EL2EDGE         = reshape(ib,6,nel); % element to bar connectivity after doubles have been merged (removed)

% Coordinates of the nodes defining each bar's end points
thBarEnds       = reshape(GCOORD_SPH4(1,EDGE2NOD'),2,nedge);
phBarEnds       = reshape(GCOORD_SPH4(2,EDGE2NOD'),2,nedge);
rBarEnds        = reshape(GCOORD_SPH4(3,EDGE2NOD'),2,nedge);

% Check if phi angles (for bar nodes) are separated an angular distance bigger than pi 
phBarEnds       = check_ang_dist_phi(GCOORD_SPH4,EDGE2NOD,phBarEnds);

% Create two new nodes on each edge
LthBar          = thBarEnds(2,:)-thBarEnds(1,:);
thBarThirds     = [thBarEnds(1,:)+(1/3)*LthBar
                   thBarEnds(1,:)+(2/3)*LthBar];
thBarThirds     = thBarThirds(:)';
LphBar          = phBarEnds(2,:)-phBarEnds(1,:);
phBarThirds     = [phBarEnds(1,:)+(1/3)*LphBar
                   phBarEnds(1,:)+(2/3)*LphBar];
phBarThirds     = phBarThirds(:)';
LrBar           = rBarEnds(2,:)-rBarEnds(1,:);
rBarThirds      = [rBarEnds(1,:)+(1/3)*LrBar
                   rBarEnds(1,:)+(2/3)*LrBar];
rBarThirds      = rBarThirds(:)';

if ~isempty(EDGE2NOD_in_cone)
    % Compute positions of EDGE2NOD_in_cone in the EDGE2NOD list
    [~,LOCB]        = ismember(EDGE2NOD_in_cone,EDGE2NOD,'rows');
    [~,LOCB2]       = ismember(fliplr(EDGE2NOD_in_cone),EDGE2NOD,'rows'); % this is just in case some bars were flipped when using find_unique_bars subfunction
    LOCB(LOCB == 0) = LOCB2(LOCB2 ~= 0); clear LOCB2;
    LOCD = [(LOCB*2-1)'; (LOCB*2)'];
    LOCD = LOCD(:)';
    % Substitute edge nodes of elements in the cone (elements within the cone boundary + elements crossing the cone boundary)
    thBarThirds(LOCD) = GCOORD_SPH_bars_els_in_cone(1,:);
    phBarThirds(LOCD) = GCOORD_SPH_bars_els_in_cone(2,:);
    rBarThirds(LOCD)  = GCOORD_SPH_bars_els_in_cone(3,:);
end

if ~isempty(EDGE2NOD_els_cross_2pi)
    % Compute positions of EDGE2NOD_els_cross_2pi in the EDGE2NOD list
    [~,LOCB]        = ismember(EDGE2NOD_els_cross_2pi,EDGE2NOD,'rows');
    [~,LOCB2]       = ismember(fliplr(EDGE2NOD_els_cross_2pi),EDGE2NOD,'rows'); % this is just in case some bars were flipped when using find_unique_bars subfunction
    LOCB(LOCB == 0) = LOCB2(LOCB2 ~= 0); clear LOCB2;
    LOCD = [(LOCB*2-1)'; (LOCB*2)'];
    LOCD = LOCD(:)';
    % Substitute edge nodes of elements in the cone (elements within the cone boundary + elements crossing the cone boundary)
    thBarThirds(LOCD) = GCOORD_SPH_bars_els_cross_2pi(1,:);
    phBarThirds(LOCD) = GCOORD_SPH_bars_els_cross_2pi(2,:);
    rBarThirds(LOCD)  = GCOORD_SPH_bars_els_cross_2pi(3,:);
end

% storage for cubic order connectivity matrix
EL2NOD         = zeros(16,nel,'uint32'); % must be 16, face and bubble nodes added later
EL2NOD(1: 4,:) = EL2NOD4;
EL2NOD(5:16,:) = nnod4 + reshape([2*EL2EDGE(:)'-1; 2*EL2EDGE(:)'],12,[]);

% storage for quadratic order (20-node) node coordinates
nnod16                       = nnod4 + 2*nedge;
GCOORD_SPH                   = zeros(3,nnod16);
GCOORD_SPH(:,1:nnod4)        = GCOORD_SPH4;
GCOORD_SPH(1,nnod4+1:nnod16) = thBarThirds;
GCOORD_SPH(2,nnod4+1:nnod16) = phBarThirds;
GCOORD_SPH(3,nnod4+1:nnod16) = rBarThirds;

GCOORD = spherical2cartesian(GCOORD_SPH); % convert to Cartesian coordinates

% =========================================================================
% Now correct local node numbering for all edge nodes
% =========================================================================
nvertx = 4;
% local coordinates of all edge nodes (nodes 5:16)
r = [2/3 1/3  0   0   0   0  1/3 2/3 2/3 1/3  0   0 ];
s = [1/3 2/3 2/3 1/3  0   0   0   0   0   0  1/3 2/3];
t = [ 0   0  1/3 2/3 2/3 1/3  0   0  1/3 2/3  0   0 ];
N = sf_dsf_tet([r;s;t],nvertx,'matrix');

if ~isempty(ELS.els_in_cone_iso)
    els_in_cone_iso = ELS.els_in_cone_iso;
    RR_X_90_CCW     = [ 1  0  0; 0  0 -1; 0  1  0 ];   % 90° counterclockwise around the X axis rotation matrix
    GCOORD_rot      = RR_X_90_CCW * GCOORD;            % rotate coordinates
    GCOORD_SPH_rot  = cartesian2spherical(GCOORD_rot); % convert to spherical coordinates
    thn_vertx       = reshape(GCOORD_SPH_rot(1,EL2NOD(1:nvertx,els_in_cone_iso)),nvertx,[]);
    thn_check       = N' * thn_vertx;
    thn_mesh        = reshape(GCOORD_SPH_rot(1,EL2NOD(5:16,els_in_cone_iso)),12,[]);
    diffth          = abs(thn_check-thn_mesh)>6e-2;
    phn_vertx       = reshape(GCOORD_SPH_rot(2,EL2NOD(1:nvertx,els_in_cone_iso)),nvertx,[]);
    phn_check       = N' * phn_vertx;
    phn_mesh        = reshape(GCOORD_SPH_rot(2,EL2NOD(5:16,els_in_cone_iso)),12,[]);
    diffph          = abs(phn_check-phn_mesh)>6e-2;
    rn_vertx        = reshape(GCOORD_SPH_rot(3,EL2NOD(1:nvertx,els_in_cone_iso)),nvertx,[]);
    rn_check        = N' * rn_vertx;
    rn_mesh         = reshape(GCOORD_SPH_rot(3,EL2NOD(5:16,els_in_cone_iso)),12,[]);
    diffr           = abs(rn_check-rn_mesh)>6e-2;
    for inod=1:2:12
        iel_flip = find(diffth(inod,:) | diffph(inod,:) | diffr(inod,:));
        if ~isempty(iel_flip)
            EL2NOD([inod+4 inod+5],els_in_cone_iso(iel_flip)) = EL2NOD([inod+5 inod+4],els_in_cone_iso(iel_flip));
        end
    end
end

if ~isempty(ELS.els_in_cone_no_iso)
    els_in_cone_no_iso = ELS.els_in_cone_no_iso;
    RR_X_90_CCW        = [ 1  0  0; 0  0 -1; 0  1  0 ];   % 90° counterclockwise around the X axis rotation matrix
    GCOORD_rot         = RR_X_90_CCW * GCOORD;            % rotate coordinates
    GCOORD_SPH_rot     = cartesian2spherical(GCOORD_rot); % convert to spherical coordinates
    thn_vertx          = reshape(GCOORD_SPH_rot(1,EL2NOD(1:nvertx,els_in_cone_no_iso)),nvertx,[]);
    thn_check          = N' * thn_vertx;
    thn_mesh           = reshape(GCOORD_SPH_rot(1,EL2NOD(5:16,els_in_cone_no_iso)),12,[]);
    diffth             = abs(thn_check-thn_mesh)>1e-12;
    phn_vertx          = reshape(GCOORD_SPH_rot(2,EL2NOD(1:nvertx,els_in_cone_no_iso)),nvertx,[]);
    phn_check          = N' * phn_vertx;
    phn_mesh           = reshape(GCOORD_SPH_rot(2,EL2NOD(5:16,els_in_cone_no_iso)),12,[]);
    diffph             = abs(phn_check-phn_mesh)>1e-12;
    rn_vertx           = reshape(GCOORD_SPH_rot(3,EL2NOD(1:nvertx,els_in_cone_no_iso)),nvertx,[]);
    rn_check           = N' * rn_vertx;
    rn_mesh            = reshape(GCOORD_SPH_rot(3,EL2NOD(5:16,els_in_cone_no_iso)),12,[]);
    diffr              = abs(rn_check-rn_mesh)>1e-08;
    for inod=1:2:12
        iel_flip = find(diffth(inod,:) | diffph(inod,:) | diffr(inod,:));
        if ~isempty(iel_flip)
            EL2NOD([inod+4 inod+5],els_in_cone_no_iso(iel_flip)) = EL2NOD([inod+5 inod+4],els_in_cone_no_iso(iel_flip));
        end
    end
end

if ~isempty(ELS.els_out_cone_cross_2pi)
    els_out_cone_cross_2pi = ELS.els_out_cone_cross_2pi;
    RR_Z_180_CCW           = [-1  0  0;   0 -1  0; 0  0  1 ]; % 180° counterclockwise around the Z axis rotation matrix
    GCOORD_rot             = RR_Z_180_CCW * GCOORD;           % rotate coordinates
    GCOORD_SPH_rot         = cartesian2spherical(GCOORD_rot); % convert to spherical coordinates
    thn_vertx              = reshape(GCOORD_SPH_rot(1,EL2NOD(1:nvertx,els_out_cone_cross_2pi)),nvertx,[]);
    thn_check              = N' * thn_vertx;
    thn_mesh               = reshape(GCOORD_SPH_rot(1,EL2NOD(5:16,els_out_cone_cross_2pi)),12,[]);
    diffth                 = abs(thn_check-thn_mesh)>1e-12;
    phn_vertx              = reshape(GCOORD_SPH_rot(2,EL2NOD(1:nvertx,els_out_cone_cross_2pi)),nvertx,[]);
    phn_check              = N' * phn_vertx;
    phn_mesh               = reshape(GCOORD_SPH_rot(2,EL2NOD(5:16,els_out_cone_cross_2pi)),12,[]);
    diffph                 = abs(phn_check-phn_mesh)>1e-12;
    rn_vertx               = reshape(GCOORD_SPH_rot(3,EL2NOD(1:nvertx,els_out_cone_cross_2pi)),nvertx,[]);
    rn_check               = N' * rn_vertx;
    rn_mesh                = reshape(GCOORD_SPH_rot(3,EL2NOD(5:16,els_out_cone_cross_2pi)),12,[]);
    diffr                  = abs(rn_check-rn_mesh)>1e-08;
    for inod=1:2:12
        iel_flip = find(diffth(inod,:) | diffph(inod,:) | diffr(inod,:));
        if ~isempty(iel_flip)
            EL2NOD([inod+4 inod+5],els_out_cone_cross_2pi(iel_flip)) = EL2NOD([inod+5 inod+4],els_out_cone_cross_2pi(iel_flip));
        end
    end
end

if ~isempty(ELS.els_out_cone_no_cross_2pi)
    els_out_cone_no_cross_2pi = ELS.els_out_cone_no_cross_2pi;
    thn_vertx                 = reshape(GCOORD_SPH(1,EL2NOD(1:nvertx,els_out_cone_no_cross_2pi)),nvertx,[]);
    thn_check                 = N' * thn_vertx;
    thn_mesh                  = reshape(GCOORD_SPH(1,EL2NOD(5:16,els_out_cone_no_cross_2pi)),12,[]);
    diffth                    = abs(thn_check-thn_mesh)>1e-12;
    phn_vertx                 = reshape(GCOORD_SPH(2,EL2NOD(1:nvertx,els_out_cone_no_cross_2pi)),nvertx,[]);
    phn_check                 = N' * phn_vertx;
    phn_mesh                  = reshape(GCOORD_SPH(2,EL2NOD(5:16,els_out_cone_no_cross_2pi)),12,[]);
    diffph                    = abs(phn_check-phn_mesh)>1e-12;
    rn_vertx                  = reshape(GCOORD_SPH(3,EL2NOD(1:nvertx,els_out_cone_no_cross_2pi)),nvertx,[]);
    rn_check                  = N' * rn_vertx;
    rn_mesh                   = reshape(GCOORD_SPH(3,EL2NOD(5:16,els_out_cone_no_cross_2pi)),12,[]);
    diffr                     = abs(rn_check-rn_mesh)>1e-12;
    for inod=1:2:12
        iel_flip = find(diffth(inod,:) | diffph(inod,:) | diffr(inod,:));
        if ~isempty(iel_flip)
            EL2NOD([inod+4 inod+5],els_out_cone_no_cross_2pi(iel_flip)) = EL2NOD([inod+5 inod+4],els_out_cone_no_cross_2pi(iel_flip));
        end
    end
end

%=========================================================================================================================
% COMPUTE A NODE ON EACH OF THE FACES OF THE TETRAHEDRA
%=========================================================================================================================
[~,EL2NOD,FACE2NOD]   = tetmesh_add_facenods(GCOORD_SPH,EL2NOD);
nface                 = size(FACE2NOD,1);
th_face_nodes_all_els = sum(reshape(GCOORD_SPH(1,FACE2NOD'),3,nface))./3;
ph_face_nodes_all_els = sum(reshape(GCOORD_SPH(2,FACE2NOD'),3,nface))./3;
r_face_nodes_all_els  = sum(reshape(GCOORD_SPH(3,FACE2NOD'),3,nface))./3;

if ~isempty(FACE2NOD_in_cone)
    % Compute positions of FACE2NOD_in_cone in the FACE2NOD list taking into account the possible permutations
    face_perms               = perms([3 2 1]);
    [~,LOCF]                 = ismember(FACE2NOD_in_cone(:,face_perms(1,:)),FACE2NOD,'rows');
    for i = 2:6
        [~,LOCF2]            = ismember(FACE2NOD_in_cone(:,face_perms(i,:)),FACE2NOD,'rows');
        LOCF(LOCF2 ~= 0)     = LOCF2(LOCF2 ~= 0);
    end
    th_face_nodes_all_els(LOCF) = GCOORD_SPH_face_nodes_els_in_cone(1,:);
    ph_face_nodes_all_els(LOCF) = GCOORD_SPH_face_nodes_els_in_cone(2,:);
    r_face_nodes_all_els(LOCF)  = GCOORD_SPH_face_nodes_els_in_cone(3,:);
end

if ~isempty(FACE2NOD_els_cross_2pi)
    % Compute positions of FACE2NOD_els_cross_2pi in the FACE2NOD list taking into account the possible permutations
    face_perms               = perms([3 2 1]);
    [~,LOCF]                 = ismember(FACE2NOD_els_cross_2pi(:,face_perms(1,:)),FACE2NOD,'rows');
    for i = 2:6
        [~,LOCF2]            = ismember(FACE2NOD_els_cross_2pi(:,face_perms(i,:)),FACE2NOD,'rows');
        LOCF(LOCF2 ~= 0)     = LOCF2(LOCF2 ~= 0);
    end
    th_face_nodes_all_els(LOCF) = GCOORD_SPH_face_nodes_els_cross_2pi(1,:);
    ph_face_nodes_all_els(LOCF) = GCOORD_SPH_face_nodes_els_cross_2pi(2,:);
    r_face_nodes_all_els(LOCF)  = GCOORD_SPH_face_nodes_els_cross_2pi(3,:);
end

GCOORD_SPH = [GCOORD_SPH [th_face_nodes_all_els; ph_face_nodes_all_els; r_face_nodes_all_els]];

GCOORD     = spherical2cartesian(GCOORD_SPH); % convert to Cartesian coordinates

end % END OF SUBFUNCTION compute_cubic_nodes_sph_rot_frame_all_els

% #####################################################################################################

function plot_points_and_elements(GCOORD,GCOORD_SPH,EL2NOD,gX_PT,gTH_PT,els,theta_cone,color_els,color_nodes,Fig_num)

figure(Fig_num)
clf
subplot(1,2,1)
% Plot the cone around z axis (point at the origin)
lightMag  = 0.80*[1 0.7 1]; % colour for the cone boundary
lightGrey = 0.90*[1 1 1]; % colour for the shell boundaries
th        = theta_cone*pi/180;
h         = 6371;
r         = h*tan(th);
[R,A]     = meshgrid(linspace(0,r,2),linspace(0,2*pi,21));
X         = R .* cos(A);
Y         = R .* sin(A);
Z         = R/tan(th);
surface(X,Y,Z,'FaceColor','none','EdgeColor',lightMag)
hold on
surface(X,Y,-Z,'FaceColor','none','EdgeColor',lightMag)
% Plot the shell
[x_sph,y_sph,z_sph] = sphere(20);
x_sph               = x_sph*3471;
y_sph               = y_sph*3471;
z_sph               = z_sph*3471;
axis equal
surface(x_sph,y_sph,z_sph,'FaceColor','none','EdgeColor',lightGrey)
[x_sph,y_sph,z_sph] = sphere(30);
x_sph               = x_sph*6371;
y_sph               = y_sph*6371;
z_sph               = z_sph*6371;
surface(x_sph,y_sph,z_sph,'FaceColor','none','EdgeColor',lightGrey)
axis([-6371 6371 -6371 6371 -6371 6371])
view(142.5,30)
grid on
hold on
tetramesh(EL2NOD(1:4,els)',GCOORD','FaceColor',color_els,'FaceAlpha',0.3)
scatter3(gX_PT(1,:),gX_PT(2,:),gX_PT(3,:),'MarkerEdgeColor','k','MarkerFaceColor',color_nodes)
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')

subplot(1,2,2)
GCOORD_SPH_plot      = GCOORD_SPH;
GCOORD_SPH_plot(3,:) = GCOORD_SPH_plot(3,:)/1000;
hold on
axis equal
axis([0 pi 0 2*pi 3 6.5])
view(142.5,30)
grid on
tetramesh(EL2NOD(1:4,els)',GCOORD_SPH_plot','FaceColor',color_els,'FaceAlpha',0.3)
scatter3(gTH_PT(1,:),gTH_PT(2,:),gTH_PT(3,:)/1e3,'MarkerEdgeColor','k','MarkerFaceColor',color_nodes)
xlabel('\theta (rad)')
ylabel('\phi (rad)')
zlabel('r (10^3 km)')

end % END OF SUBFUNCTION plot_points_and_elements

% #####################################################################################################

function PointID = calc_PointIDs(GCOORD,EL2NOD,PointID_c,DB_indices,FigNo)

nVnod = max(max(EL2NOD(1: 4,:)));
nnod  = max(max(EL2NOD(5:10,:)));
nEnod = nnod-nVnod;

% Convert to old format (have to rewrite this function to use
% PointID-format)
% ===========================================================
DBnods_c      = find(PointID_c>0)';
if isempty(DBnods_c)
    PointID   = zeros(1,nnod,'int32');
    return
end
DBnods_c(:,2) = PointID_c( DBnods_c );

% Get boundary index information for new nodes
% ================================================
% 1) Create a pointer that for each edge node in the mesh returns the 2
%    bounding vertex nodes
pointer_Edge2Vnod                       = zeros(nEnod,2,'int32');
pointer_Edge2Vnod(EL2NOD( 5,:)-nVnod,:) = EL2NOD([1 2],:)';
pointer_Edge2Vnod(EL2NOD( 6,:)-nVnod,:) = EL2NOD([2 3],:)';
pointer_Edge2Vnod(EL2NOD( 7,:)-nVnod,:) = EL2NOD([3 4],:)';
pointer_Edge2Vnod(EL2NOD( 8,:)-nVnod,:) = EL2NOD([4 1],:)';
pointer_Edge2Vnod(EL2NOD( 9,:)-nVnod,:) = EL2NOD([1 3],:)';
pointer_Edge2Vnod(EL2NOD(10,:)-nVnod,:) = EL2NOD([2 4],:)';


% (2) Find fine mesh vertex nodes (i.e. coarse mesh nodes) that are in the
%     list of coarse domain boundary nodes (1st column in DBnods_c)
[tf,loc] = ismember(pointer_Edge2Vnod,DBnods_c(:,1)); % find DBnod_c in pointer_Edge2Vnod
tf       = find(tf);
% Create new pointer to the domain boundary index (2nd column in DBnods_c)
pointer_Edge2DBnod     = zeros(nEnod,2,'int32'); % create new pointer
pointer_Edge2DBnod(tf) = DBnods_c(loc(tf),2);

% (3)  Use the pointer pointer_Edge2DBnod to check (for each NEW node) 
%      what the domain boundary index of its two neighborung vertex nodes
%      are. The domain boundary index for the new nodes can be found using
%      the following logic:
%      Edge nodes cannot be generated in domain corners but only on
%      domain edges or faces.
%      Face nodes cannot be generated on domain edges (and of course not in
%      domain corners).
%
%    108---211---107    % Coordinate system to define xmin, xmax, ymin, ymax, zmin, zmax
%     /:         /|     %     y
%  212 :  306  210|     %   /
%   /  :       /  |     %  /
%105---:209---106 |     % 0------- x
%  |  208  304|  207    % |
%  |   :      |   |     % |
%  |305:      |303|     % |
% 205  : 302 206  |     %-z
%  | 104---203|---103
%  |  /       |  /
%  |204  301  |202
%  |/         |/
%101--- 201---102

% Faces: 301 bottom
%        302 front
%        303 right
%        304 back
%        305 left
%        306 top
% Edges: 201 lower front
%        202 lower right
%        203 lower back
%        204 lower left
%        205 front left
%        206 front right
%        207 back right
%        208 back left
%        209 top front
%        210 top right
%        211 top back
%        212 top left
% Corners: 101 lower front left
%          102 lower front right
%          103 lower back right
%          104 lower back left
%          105 top front left
%          106 top front right
%          107 top back right
%          108 top back left
%
% special: 310 face nodes west of MAR
%          220 edge nodes west of MAR
%          311 face nodes east of MAR
%          221 edge nodes east of MAR
%          222 nodes on MAR
%          223 nodes in TF
%          110 nodes at JCT (intersection MAR-TF)

% There are 9 standard cases and 6 special cases in 3D:
%       1) both on same F ==> F
%       2) both on same E ==> E
%       3) C and E ==> E (C+E must be on same domain edge)
%       4) C and F ==> F (C+F must be on same domain face)
%       5) E and F ==> F (E+F must be on same domain face)
%       6) C and 0 ==> 0
%       7) E and 0 ==> 0
%       8) F and 0 ==> 0
%       9) 0 and 0 ==> 0 
%      S1) C1 and C2 ==> E between C1 and C2
%      S2) E1 and E2 ==> F between E1 and E2
%      S3) F1 and F2 ==> 0 (inside domain)
%      S4) C1 and E2 ==> 0 (inside domain)
%      S5) C1 and F2 ==> 0 (inside domain)
%      S6) E1 and F2 ==> 0 (inside domain)
%      The numbers in the special cases indicate DIFFERENT domain faces !!!

% % % For testing and understanding the logic:
% % pointer_Edge2DBnod = [301 301; % case 1
% %                        201 201; % case 2
% %                        101 201; % case 3
% %                        101 301; % case 4
% %                        201 301; % case 5
% %                        101  0 ; % case 6
% %                        201  0 ; % case 7
% %                        301  0 ; % case 8
% %                         0   0 ; % case 9
% %                        101 102; % special case 1
% %                        201 202; % special case 2
% %                        301 302];% special case 3
% %                        101 202; % special case 4
% %                        101 303; % special case 5
% %                        201 303];% special case 6

% Cases 1:9 are taken care of by chosing the LARGER index of the 2
% bounding vertex nodes. Treat all nodes as standard cases first:
DBindx  = max( pointer_Edge2DBnod,[],2 );
% HOWEVER, if any index is zero, the node will be inside the domain. Set
% these indices equal to zero now.
ind_0   = any(pointer_Edge2DBnod==0,2);
DBindx(ind_0) = 0;

% Need special treatment for the special cases...
% (a) all potential special cases
indS  = find( pointer_Edge2DBnod(:,1)~=pointer_Edge2DBnod(:,2) & ...
              pointer_Edge2DBnod(:,1)>0 & pointer_Edge2DBnod(:,2)>0 );
% (b) special case 1 (between 2 corners); all possible combinations
indS1 = indS(all(pointer_Edge2DBnod(indS,:)>100,2) & all(pointer_Edge2DBnod(indS,:)<200,2));
if ~isempty(indS1)
    ind_0          = indS1;
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[101 102]) ,2) );
    DBindx(indS1a) = 201; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[102 103]) ,2) );
    DBindx(indS1a) = 202; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[103 104]) ,2) );
    DBindx(indS1a) = 203; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[101 104]) ,2) );
    DBindx(indS1a) = 204; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[101 105]) ,2) );
    DBindx(indS1a) = 205; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[102 106]) ,2) );
    DBindx(indS1a) = 206; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[103 107]) ,2) );
    DBindx(indS1a) = 207; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[104 108]) ,2) );
    DBindx(indS1a) = 208; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[105 106]) ,2) );
    DBindx(indS1a) = 209; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[106 107]) ,2) );
    DBindx(indS1a) = 210; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[107 108]) ,2) );
    DBindx(indS1a) = 211; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[105 108]) ,2) );
    DBindx(indS1a) = 212; ind_0 = setdiff(ind_0,indS1a);
    
    for iDface=1:length(DB_indices)
        iSameFace = sum(ismember( pointer_Edge2DBnod(ind_0,:),DB_indices{iDface} ),2)==2;
        DBindx(ind_0(iSameFace)) = max(DB_indices{iDface});
        ind_0(iSameFace)         = [];
    end
    
    DBindx(ind_0)  = 0;
end

% (c) special case 2 (between 2 different edges); all possible combinations
indS2 = indS(all(pointer_Edge2DBnod(indS,:)>200,2) & all(pointer_Edge2DBnod(indS,:)<300,2));
if ~isempty(indS2)
    ind_0          = indS2;
    indS2a         = indS2( all( ismember(pointer_Edge2DBnod(indS2,:),[201 202 203 204]) ,2) );
    DBindx(indS2a) = 301; ind_0 = setdiff(ind_0,indS2a);
    indS2a         = indS2( all( ismember(pointer_Edge2DBnod(indS2,:),[201 205 206 209]) ,2) );
    DBindx(indS2a) = 302; ind_0 = setdiff(ind_0,indS2a);
    indS2a         = indS2( all( ismember(pointer_Edge2DBnod(indS2,:),[202 206 207 210]) ,2) );
    DBindx(indS2a) = 303; ind_0 = setdiff(ind_0,indS2a);
    indS2a         = indS2( all( ismember(pointer_Edge2DBnod(indS2,:),[203 207 208 211]) ,2) );
    DBindx(indS2a) = 304; ind_0 = setdiff(ind_0,indS2a);
    indS2a         = indS2( all( ismember(pointer_Edge2DBnod(indS2,:),[204 205 208 212]) ,2) );
    DBindx(indS2a) = 305; ind_0 = setdiff(ind_0,indS2a);
    indS2a         = indS2( all( ismember(pointer_Edge2DBnod(indS2,:),[209 210 211 212]) ,2) );
    DBindx(indS2a) = 306; ind_0 = setdiff(ind_0,indS2a);
    
    for iDface=1:length(DB_indices)
        iSameFace = sum(ismember( pointer_Edge2DBnod(ind_0,:),DB_indices{iDface} ),2)==2;
        DBindx(ind_0(iSameFace)) = max(DB_indices{iDface});
        ind_0(iSameFace)         = [];
    end
    
    DBindx(ind_0)  = 0;
end

% (d) special case 3
indS3 = indS(all(pointer_Edge2DBnod(indS,:)>300,2) & all(pointer_Edge2DBnod(indS,:)<400,2));

% (e) special cases 4, 5 and 6
indS456 = setdiff(indS,[indS1(:);indS2(:);indS3(:)]);
% Now loop over the domain faces and check if the different DB indices
% (e.g. C1 and E2) are on the same domain face as defined by the cell-matrix 
% "DB_indices". If they are on the same face, they are NOT INSIDE the
% domain!
for iDface=1:length(DB_indices)
    iSameFace = sum(ismember( pointer_Edge2DBnod(indS456,:),DB_indices{iDface} ),2)==2;
    indS456(iSameFace) = [];
end
% old version:
% Dfaces  = [301 201 202 203 204 101 102 103 104; % domain bottom     (zmin)
%            302 201 205 206 209 101 102 105 106; % domain front      (ymin)
%            303 202 206 207 210 102 103 106 107; % domain right side (xmax)
%            304 203 207 208 211 103 104 107 108; % domain back       (ymax)
%            305 204 205 208 212 101 104 105 108; % domain left side  (xmin)
%            306 209 210 211 212 105 106 107 108];% domain top        (zmax)
% % Now loop over the domain faces and check if the different DB indices
% % (e.g. C1 and E2) are on the same domain face as defined by matrix Dfaces.
% for iDface=1:size(Dfaces,1)
%     iSameFace = sum(ismember( pointer_Edge2DBnod(indS456,:),Dfaces(iDface,:) ),2)==2;
%     indS456(iSameFace) = [];
% end

% Special case 3,4,5 and 6 (between 2 different domain faces thus INSIDE
% the domain) ===> simply remove these nodes from the list of boundary nodes
inod_DB = (nVnod+1:nnod)';
DBindx ([indS3(:);indS456(:)]) = [];
inod_DB([indS3(:);indS456(:)]) = [];

% Remove all nodes that have a zero DB index
ind_0          = ~DBindx;
DBindx(ind_0)  = [];
inod_DB(ind_0) = [];

% merge coarse-mesh and new domain boundary nodes
DBnods                = [DBnods_c; [inod_DB DBindx]];
PointID              = zeros(1,nnod,'int32');
PointID(DBnods(:,1)) = DBnods(:,2);

% FigNo = 11;
if FigNo
    sfigure(FigNo);clf;
    
    uniqDBindx        = unique(DBnods(:,2));
    nDBindx           = length(uniqDBindx);
    col               = zeros(max(DBnods(:,2)),3);
    col(uniqDBindx,:) = lines(nDBindx);

    tetramesh(EL2NOD(1:4,:)',GCOORD',...
              'EdgeColor','g',...
              'FaceAlpha',0.6,'FaceColor','k',...
              'LineStyle','-','LineWidth',1);
    hold all
    xlabel('X');ylabel('Y');zlabel('Z');
    for ii=1:length(DBnods)
        nod  = DBnods(ii,1);
        indx = DBnods(ii,2);
        text(GCOORD(1,nod),GCOORD(2,nod),GCOORD(3,nod),num2str(indx),...
             'Fontweight','bold','Fontsize',14,'Color',col(indx,:));
        scatter3(GCOORD(1,nod),GCOORD(2,nod),GCOORD(3,nod),5,...
                 'Marker','o','LineWidth',2,...
                 'MarkerFaceColor',col(indx,:),'MarkerEdgeColor',col(indx,:));
    end
    
    sfigure(FigNo+1);clf;
    tetramesh(EL2NOD(1:4,:)',GCOORD',...
              'EdgeColor','g',...
              'FaceAlpha',0.6,'FaceColor','k',...
              'LineStyle','-','LineWidth',1);
    hold all
    xlabel('X');ylabel('Y');zlabel('Z');
    for ii=1:length(DBnods)
        nod  = DBnods(ii,1);
        indx = DBnods(ii,2);
        text(GCOORD(1,nod),GCOORD(2,nod),GCOORD(3,nod),num2str(nod),...
             'Fontweight','bold','Fontsize',14,'Color',col(indx,:));
        scatter3(GCOORD(1,nod),GCOORD(2,nod),GCOORD(3,nod),5,...
                 'Marker','o','LineWidth',2,...
                 'MarkerFaceColor',col(indx,:),'MarkerEdgeColor',col(indx,:));
    end
end

end % END OF SUBFUNCTION calc_PointIDs