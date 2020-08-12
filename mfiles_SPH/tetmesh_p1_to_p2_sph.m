function [GCOORD_curved,GCOORD_SPH,EL2NOD,PointID,els_out_cone_no_cross_2pi,els_out_cone_cross_2pi,els_in_cone_no_iso,els_in_cone_iso] = ...
    tetmesh_p1_to_p2_sph(GCOORD_c,GCOORD_SPH_c,EL2NOD_c,PointID_c,DB_indices,r_cmb,r_surf,SETTINGS)
% Usage: [GCOORD_curved,GCOORD_SPH,EL2NOD,PointID,els_out_cone_no_cross_2pi,els_out_cone_cross_2pi,els_in_cone_no_iso,els_in_cone_iso] = ...
%   tetmesh_p1_to_p2_sph(GCOORD_c,GCOORD_SPH_c,EL2NOD_c,PointID_c,DB_indices,r_cmb,r_surf,SETTINGS)
%
% Purpose: 
%   Creates a quadratic order (10-node) from a linear (4-node) tetrahedra
%   mesh by adding edge nodes to all elements. 
%   To do so, bars are created that connect the vertex nodes between which
%   a new node has to be created. Then a unique list of these bars is
%   calculated to avoid multiple nodes in the same place. Once the unique
%   bar-list is obtained, node coordinates and the new connectivity matrix
%   are calculated.
%   
%   This routine returns 10-node tetrahedrons with curved edges since the
%   edges are being splitted in spherical coordinates. Remember that a
%   straight edge in spherical coordinates (theta,phi,r) means a curved
%   edge in Cartesian coordinates (x,y,z).
%   However, elements inside a cone (given by an angle ~45°) are very
%   distorted when mid-side nodes are computed in spherical coordinates. 
%   This is mainly due to two facts:
%       1) Elements that are crossed by Z axis. These elements present a 
%       wide longitude (phi) range for the vertices within the same element
%       leading to inaccurate mid-side positions.
%       2) Elements with 1 or 2 nodes on the Z axis. Any node on the Z axis
%       does not have defined the longitude (phi). However MATLAB assigns 0
%       value to the longitude of those nodes.
%   Elements crossing phi = 2pi also show issues when computing mid-side
%   nodes. Their vertices might be close in Cartesian coordinates, however, 
%   they are ~2pi radians (360 degrees) far from each other in spherical
%   coordinates. 
%   
%   In order to deal with these two issues it is needed a classification of
%   elements in function on their position in the mesh: 
%   - els_out_cone_no_cross_2pi:
%       Elements outside the cone ( = 4 vertices outside the cone) and not
%       crossing phi = 2pi. These elements have straight edges in the
%       original spherical frame. 
%   - els_out_cone_cross_2pi:
%       Elements outside the cone ( = 4 vertices outside the cone) and
%       crossing phi = 2pi. These elements have straight edges in the
%       rotated spherical frame (180° around Z axis). 
%   - els_in_cone_no_iso:
%       Elements within the cone ( = 4 vertices within the cone) + elements
%       crossing the cone boundary and having only 1 vertex outside the
%       cone ( = 3 vertices within the cone + 1 vertex outside the cone).
%       These elements have straight edges in the rotated spherical frame
%       (90° around X axis). 
%   - els_in_cone_iso:
%       Elements crossing the cone boundary and having at least 1 edge
%       (bar) outside the cone ( = 2 vertices within the cone + 2 vertices
%       outside the cone or 1 vertex within the cone + 3 vertices outside
%       the cone). These elements have curved edges in the rotated
%       spherical frame (90° around X axis) since: 
%           (1) some of their mid-side nodes are computed in the original
%               spherical frame (bar or bars that are outside the cone ->
%               straight edges in the original spherical frame).
%           (2) some others mid-side nodes are computed in the rotated
%               spherical frame (bar or bars that are inside the cone or
%               crossing the cone boundary -> straight edges in the rotated
%               spherical frame (90° around X axis)). 
%       Note that the number of elements crossing the cone boundary DOES 
%       NOT HAVE TO BE THE SAME as the number of elements crossing the cone
%       boundary and having at least 1 edge (bar) outside the cone.
%   
%   Steps to compute mid-side nodes:
%   1. Compute mid-side nodes for elements outside the cone and crossing
%      phi = 2pi in a rotated spherical frame (180° around Z axis):
%       - Make a 180° counterclockwise rotation around the Z axis for the
%         elements outside the cone and crossing phi = 2pi
%         (els_out_cone_cross_2pi).
%       - Convert to spherical coordinates.
%       - Compute bars for those elements (bar2node_els_out_cone_cross_2pi)
%       - Compute mid-side nodes in straight edges of the rotated spherical
%         frame.
%       - Convert back to Cartesian coordinates.
%       - Undo the rotation (make a 180° clockwise rotation around the Z
%         axis).
%       - Convert again to spherical coordinates
%         (GCOORD_SPH_els_out_cone_cross_2pi).
%   
%   2. Compute mid-side nodes for elements in the cone (els_in_cone) in a 
%      rotated spherical frame (90° around X axis): 
%       - Make a 90° counterclockwise rotation around the X axis for the
%         elements in the cone (els_in_cone = els_in_cone_no_iso + 
%         els_in_cone_iso). 
%       - Convert to spherical coordinates.
%       - Compute bars for elements in the cone (bar2node_els_in_cone).
%       - Find bars for those elements that have at least 1 edge (bar)
%         outside the cone
%         (bar2node_els_in_cone_isoparametric_with_2_bars_outside_cone).
%       - Remove these outside cone bars in the bar2node_els_in_cone list.
%       - Compute mid-side nodes in straight edges of the rotated spherical
%         frame.
%       - Convert back to Cartesian coordinates.
%       - Undo the rotation (make a 90° clockwise rotation around the X
%         axis).
%       - Convert again to spherical coordinates (GCOORD_SPH_els_in_cone).
%   
%   3. Compute mid-side nodes for all elements in the original spherical
%      frame: 
%       - Compute bars for all elements (bar2node).
%       - Compute mid-side nodes in straight edges of the original
%         spherical frame.
%       - Find bar positions of bar2node_els_in_cone in bar2node.
%       - Substitute GCOORD_SPH_els_in_cone in the mid-side nodes computed
%         in straight edges of the original spherical frame using the bar
%         positions of bar2node_els_in_cone in bar2node.
%       - Find bar positions of bar2node_els_out_cone_crossing_2pi in
%         bar2node.
%       - Substitute GCOORD_SPH_els_out_cone_crossing_2pi in the mid-side
%         nodes computed in straight edges of the original spherical frame
%         using the bar positions of bar2node_els_out_cone_crossing_2pi in
%         bar2node.
%
% Input:
%   GCOORD_c     : [matrix]    : Cartesian coordinates of all nodes 
%                                (4 nodel) (3 x nnod_c)
%   GCOORD_SPH_c : [matrix]    : spherical coordinates of all nodes 
%                                (4 nodel) (3 x nnod_c)
%   EL2NOD_c     : [matrix]    : connectivity matrix (4 x nel)
%   PointID_c    : [vector]    : domain boundary index for each node
%   DB_indices   : [matrix]    : domain boundary indices
%   r_cmb        : [scalar]    : radius for CMB (km)
%   r_surf       : [scalar]    : radius for surface (km)
%   SETTINGS     : [structure] : model parameters
%
% Output:
%   GCOORD_curved             : [matrix]   : Cartesian coordinates of all
%                                            nodes new mesh (3 x ?)
%   GCOORD_SPH                : [matrix]   : spherical coordinates of all
%                                            nodes new mesh (3 x ?)
%   EL2NOD                    : [matrix]   : connectivity matrix (10 x nel)
%   PointID                   : [vector]   : domain boundary index for each
%                                            node (1 x ?)
%   els_out_cone_no_cross_2pi : [vector]   : elements outside the cone (4 
%                                            vertices outside the cone) and
%                                            not crossing phi = 2pi. These
%                                            elements have straight edges 
%                                            in the original spherical
%                                            frame
%   els_out_cone_cross_2pi    : [vector]   : elements outside the cone (4 
%                                            vertices outside the cone) and
%                                            crossing phi = 2pi. These
%                                            elements have straight edges 
%                                            in the rotated spherical frame
%                                            (180° around Z axis)
%   els_in_cone_no_iso        : [vector]   : elements within the cone 
%                                            (4 vertices within the cone) +
%                                            elements crossing the cone
%                                            boundary and having only 1
%                                            vertex outside the cone 
%                                            (3 vertices within the cone +
%                                            1 vertex outside the cone) 
%                                            These elements have straight
%                                            edges in the rotated spherical
%                                            frame (90° around X axis). 
%   els_in_cone_iso           : [vector]   : elements crossing the cone
%                                            boundary and having at least 1
%                                            edge (bar) outside the cone 
%                                            (2 vertices within the cone + 
%                                            2 vertices outside the cone OR 
%                                            1 vertex within the cone + 
%                                            3 vertices outside the cone).
%                                            These elements have curved
%                                            edges in the rotated spherical
%                                            frame (90° around X axis)
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH March 2013
% JMT Aug 2016: Compute curved edges for a spherical mesh
%

%======================================================================================================
% COMPUTE ELEMENTS IN RELATION WITH THE CONE AND CROSSING PHI = 2PI
%======================================================================================================
[els_in_cone,els_within_cone,els_cone_bnd,els_out_cone,els_in_cone_iso] = check_els_in_cone(GCOORD_SPH_c,EL2NOD_c,SETTINGS.theta_cone);
[els_cross_2pi,~]         = check_phi(GCOORD_SPH_c,EL2NOD_c(:,els_out_cone));
els_out_cone_cross_2pi    = els_out_cone(els_cross_2pi);
els_out_cone_no_cross_2pi = els_out_cone;
els_out_cone_no_cross_2pi(ismember(els_out_cone_no_cross_2pi,els_out_cone_cross_2pi)) = [];
els_in_cone_no_iso        = els_in_cone;
els_in_cone_no_iso(ismember(els_in_cone_no_iso,els_in_cone_iso)) = [];

%======================================================================================================
% COMPUTE MID-SIDE NODES IN SPHERICAL COORDINATES FOR ELEMENTS OUTSIDE THE CONE 
% AND CROSSING PHI = 2PI AFTER A 180° ROTATION AROUND Z AXIS 
%======================================================================================================
[GCOORD_SPH_els_out_cone_cross_2pi,bar2node_els_out_cone_cross_2pi] = ...
    mid_side_nodes_sph_rot_frame_els_out_cone_cross_2pi(GCOORD_c,EL2NOD_c,els_out_cone_cross_2pi);

%======================================================================================================
% COMPUTE MID-SIDE NODES IN SPHERICAL COORDINATES FOR ELEMENTS IN THE CONE (WITHIN THE   
% CONE BOUNDARY + CROSSING THE CONE BOUNDARY) AFTER A 90° ROTATION AROUND X AXIS 
%======================================================================================================
[GCOORD_SPH_els_in_cone,bar2node_els_in_cone] = ...
    mid_side_nodes_sph_rot_frame_els_in_cone(GCOORD_c,GCOORD_SPH_c,EL2NOD_c,els_in_cone,els_in_cone_iso,SETTINGS.theta_cone);

%======================================================================================================
% COMPUTE MID-SIDE NODES IN SPHERICAL COORDINATES (CURVED EDGES) FOR ALL ELEMENTS
%======================================================================================================
[GCOORD_curved,GCOORD_SPH,EL2NOD] = ...
    mid_side_nodes_sph_frame_all_els(GCOORD_SPH_c,EL2NOD_c,GCOORD_SPH_els_in_cone,bar2node_els_in_cone,GCOORD_SPH_els_out_cone_cross_2pi,bar2node_els_out_cone_cross_2pi);

%======================================================================================================
% PLOTS (for debugging)
%======================================================================================================
SETTINGS.show_figs = 0;
if SETTINGS.show_figs
    distortion = plot_distortion_mid_side_nodes(GCOORD_c,GCOORD_SPH_c,EL2NOD_c,EL2NOD,els_in_cone_iso);
    plot_figures(GCOORD_c,GCOORD_curved,GCOORD_SPH,EL2NOD,els_in_cone_iso,els_within_cone,els_cone_bnd,els_out_cone,els_out_cone_cross_2pi,distortion,SETTINGS.theta_cone);
end

%======================================================================================================
% DATA FOR OUTPUT
%======================================================================================================
% Calculate boundary index for the new nodes
FigNo   = 0;
PointID = calc_PointIDs(GCOORD_curved,EL2NOD,PointID_c,DB_indices,FigNo);

% Make sure that boundary nodes are exactly on the boundaries (the radius of boundary nodes can shift
% slightly (1e-8) from the actual boundary because of transforming between Cartesian and spherical coord)
GCOORD_SPH(3,PointID==301) = r_cmb;
GCOORD_SPH(3,PointID==306) = r_surf;

end % END OF FUNCTION tetmesh_p1_to_p2_sph

% #####################################################################################################
%                                           SUB-FUNCTIONS
% #####################################################################################################

function [edges_uniq,IIu,JJu] = find_unique_bars(edges)

if size(edges,2)~=2
    error('expecting a nx2 matrix');
end

% find FIRST occurrence of unique rows in edges
% use sort so that vertices defining the edges are sorted increasingly
[~,iu,ju]   = unique(sort(edges,2),'rows','first');

% create index IIu to extract edges_uniq from edges
[IIu,~,JJu] = unique(iu(ju));

% extract unique edges from input
edges_uniq  = edges(IIu,:);

end % END OF SUBFUNCTION find_unique_bars

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

function [GCOORD_SPH_els_out_cone_cross_2pi,bar2node_els_out_cone_cross_2pi] = ...
    mid_side_nodes_sph_rot_frame_els_out_cone_cross_2pi(GCOORD_c,EL2NOD_c,els_out_cone_cross_2pi)
% Usage: [GCOORD_SPH_els_out_cone_crossing_2pi,bar2node_els_out_cone_crossing_2pi] = ...
%   mid_side_nodes_sph_rot_frame_els_out_cone_cross_2pi(GCOORD_c,EL2NOD_c,els_out_cone_cross_2pi)
%
% Purpose: 
%   Compute mid-side nodes for elements outside the cone and crossing 
%   phi = 2pi in a rotated spherical frame (180° around Z axis) following
%   the next steps: 
%       - Make a 180° counterclockwise rotation around the Z axis for the
%         elements outside the cone and crossing phi = 2pi
%         (els_out_cone_cross_2pi)
%       - Convert to spherical coordinates
%       - Compute bars for those elements (bar2node_els_out_cone_cross_2pi)
%       - Compute mid-side nodes in straight edges of the rotated spherical
%         frame.
%       - Convert back to Cartesian coordinates
%       - Undo the rotation (make a 180° clockwise rotation around the Z
%         axis)
%       - Convert again to spherical coordinates
%         (GCOORD_SPH_els_out_cone_cross_2pi)
%
% Input:
%   GCOORD_c               : [matrix] : Cartesian coordinates of all nodes
%                                       (4 nodel) (3 x nnod_c)
%   EL2NOD_c               : [matrix] : connectivity matrix (4 x nel)   
%   els_out_cone_cross_2pi : [vector] : indices for those elements outside
%                                       the cone and crossing phi = 2pi 
%
% Output:
%   GCOORD_SPH_els_out_cone_cross_2pi : [matrix] : spherical coordinates of
%                                                  nodes belonging those
%                                                  elements outside the
%                                                  cone and crossing 
%                                                  phi = 2pi
%   bar2node_els_out_cone_cross_2pi   : [matrix] : bar list for elements
%                                                  outside the cone and
%                                                  crossing phi = 2pi
%
% JMT Aug 2016

% Make a 180° counterclockwise rotation around the Z axis of the elements in the cone
RR_Z_180_CCW                          = [-1  0  0 ; ...
                                          0 -1  0 ; ...
                                          0  0  1 ];
GCOORD_c_rot                          = RR_Z_180_CCW * GCOORD_c;

% Convert to spherical coordinates
GCOORD_SPH_c_rot                      = cartesian2spherical(GCOORD_c_rot);

% Compute mid-side nodes in straight edges of the spherical system
EL2NOD_c_els_out_cone_cross_2pi       = EL2NOD_c(:,els_out_cone_cross_2pi);
bar2node_els_out_cone_cross_2pi       = reshape( EL2NOD_c_els_out_cone_cross_2pi([ 1  2 ... % 1st edge of parent
                                                                                   2  3 ... % 2nd edge of parent
                                                                                   3  4 ... % 3rd edge of parent
                                                                                   4  1 ... % 4th edge of parent
                                                                                   1  3 ... % 5th edge of parent
                                                                                   4  2 ... % 6th edge of parent
                                                                                       ]',:),2,[] )';

% Find the bars that are shared by neighboring elements and return a unique list of bars.
[bar2node_els_out_cone_cross_2pi,~,~] = find_unique_bars(bar2node_els_out_cone_cross_2pi); % *SUBFUNCTION*

% Coordinates of the nodes defining each bar's end points
thetaBarEnds_els_out_cone_cross_2pi   = reshape(GCOORD_SPH_c_rot(1,bar2node_els_out_cone_cross_2pi'),2,[]);
phiBarEnds_els_out_cone_cross_2pi     = reshape(GCOORD_SPH_c_rot(2,bar2node_els_out_cone_cross_2pi'),2,[]);
rBarEnds_els_out_cone_cross_2pi       = reshape(GCOORD_SPH_c_rot(3,bar2node_els_out_cone_cross_2pi'),2,[]);

% Create new node at each bar mid-point
thetaBarMids_els_out_cone_cross_2pi   = 0.5*sum(thetaBarEnds_els_out_cone_cross_2pi,1);
phiBarMids_els_out_cone_cross_2pi     = 0.5*sum(phiBarEnds_els_out_cone_cross_2pi,1);
rBarMids_els_out_cone_cross_2pi       = 0.5*sum(rBarEnds_els_out_cone_cross_2pi,1);
GCOORD_SPH_els_out_cone_cross_2pi_rot = [thetaBarMids_els_out_cone_cross_2pi; ...
                                         phiBarMids_els_out_cone_cross_2pi;   ...
                                         rBarMids_els_out_cone_cross_2pi];
 
% Convert to Cartesian coordinates
GCOORD_els_out_cone_cross_2pi_rot     = spherical2cartesian(GCOORD_SPH_els_out_cone_cross_2pi_rot);

% Undo the rotation (rotate 180° around Z axis clockwise)
RR_Z_180_CW                           = [-1  0  0 ; ...
                                          0 -1  0 ; ...
                                          0  0  1 ];
GCOORD_els_out_cone_cross_2pi         = RR_Z_180_CW * GCOORD_els_out_cone_cross_2pi_rot;

% Convert to spherical coodinates
GCOORD_SPH_els_out_cone_cross_2pi     = cartesian2spherical(GCOORD_els_out_cone_cross_2pi);

end % END OF SUBFUNCTION mid_side_nodes_sph_rot_frame_els_out_cone_cross_2pi

% #####################################################################################################

function [GCOORD_SPH_els_in_cone,bar2node_els_in_cone] = mid_side_nodes_sph_rot_frame_els_in_cone...
    (GCOORD_c,GCOORD_SPH_c,EL2NOD_c,els_in_cone,els_in_cone_iso,theta_cone)
% Usage: [GCOORD_SPH_els_in_cone,bar2node_els_in_cone] = mid_side_nodes_sph_rot_frame_els_in_cone...
%   (GCOORD_c,GCOORD_SPH_c,EL2NOD_c,els_in_cone,els_in_cone_iso,theta_cone)
%
% Purpose: 
%   Compute mid-side nodes for elements in the cone (els_in_cone) in a
%   rotated spherical frame (90° around X axis) following the next steps:
%       - Make a 90° counterclockwise rotation around the X axis for the
%         elements in the cone (els_in_cone = els_in_cone_not_isoparametric
%         + els_in_cone_isoparametric).
%       - Convert to spherical coordinates.
%       - Compute bars for elements in the cone (bar2node_els_in_cone).
%       - Find bars for those elements that have at least 1 edge (bar)
%         outside the cone
%         (bar2node_els_in_cone_isoparametric_with_2_bars_outside_cone).
%       - Remove these outside cone bars in the bar2node_els_in_cone list.
%       - Compute mid-side nodes in straight edges of the rotated spherical
%         frame.
%       - Convert back to Cartesian coordinates.
%       - Undo the rotation (make a 90° clockwise rotation around the X
%         axis).
%       - Convert again to spherical coordinates (GCOORD_SPH_els_in_cone).
%
% Input:
%   GCOORD_c               : [matrix] : Cartesian coordinates of all nodes
%                                       (4 nodel) (3 x nnod_c)
%   GCOORD_SPH_c           : [matrix] : spherical coordinates of all nodes
%                                       (4 nodel) (3 x nnod_c)
%   EL2NOD_c               : [matrix] : connectivity matrix (4 x nel)
%   els_in_cone            : [vector] : indices for those elements in the
%                                       cone (els_in_cone =
%                                       els_in_cone_no_iso+els_in_cone_iso) 
%   els_in_cone_iso        : [vector] : indices for those elements crossing
%                                       the cone boundary and having at
%                                       least 1 edge (bar) outside the cone
%                                       boundary 
%   theta_cone             : [scalar] : half aperture of the cone 
%                                       (0° < theta_cone < 90°)
%
% Output:
%   GCOORD_SPH_els_in_cone : [matrix] : spherical coordinates of nodes
%                                       belonging those elements in the
%                                       cone (els_in_cone =
%                                       els_in_cone_not_isoparametric +
%                                       els_in_cone_isoparametric)
%   bar2node_els_in_cone   : [matrix] : bar list for elements in the cone
%                                       (els_in_cone =
%                                       els_in_cone_not_isoparametric +
%                                       els_in_cone_isoparametric)
%
% JMT Aug 2016

% Make a 90° counterclockwise rotation around the X axis of the elements in the cone
RR_X_90_CCW                = [ 1  0  0 ; ...
                               0  0 -1 ; ...
                               0  1  0 ];
GCOORD_c_rot               = RR_X_90_CCW * GCOORD_c;

% Convert to spherical coordinates
GCOORD_SPH_c_rot           = cartesian2spherical(GCOORD_c_rot);

% Compute mid-side nodes in straight edges of the spherical system
EL2NOD_c_in_cone           = EL2NOD_c(:,els_in_cone);
bar2node_els_in_cone       = reshape( EL2NOD_c_in_cone([ 1  2 ... % 1st edge of parent
                                                         2  3 ... % 2nd edge of parent
                                                         3  4 ... % 3rd edge of parent
                                                         4  1 ... % 4th edge of parent
                                                         1  3 ... % 5th edge of parent
                                                         4  2 ... % 6th edge of parent
                                                              ]',:),2,[] )';

% Find the bars that are shared by neighboring elements and return a unique list of bars.
[bar2node_els_in_cone,~,~] = find_unique_bars(bar2node_els_in_cone); % *SUBFUNCTION*

% Take the elements crossing the cone boundary and having at least 1 edge (bar) outside the cone boundary (els_in_cone_iso). 
% We need to track these elements since they will need an isoparametric spherical-to-local routine (they have curved spherical sides in the spherical rotated frame). 
% Compute the bars for those elements having at least 1 edge (bar) outside the cone boundary
EL2NOD_c_els_in_cone_iso   = EL2NOD_c(:,els_in_cone_iso);
bar2node_els_in_cone_iso   = reshape( EL2NOD_c_els_in_cone_iso([ 1  2 ... % 1st edge of parent
                                                                 2  3 ... % 2nd edge of parent
                                                                 3  4 ... % 3rd edge of parent
                                                                 4  1 ... % 4th edge of parent
                                                                 1  3 ... % 5th edge of parent
                                                                 4  2 ... % 6th edge of parent
                                                                      ]',:),2,[] )';
[bar2node_els_in_cone_iso,~,~] = find_unique_bars(bar2node_els_in_cone_iso); % *SUBFUNCTION*

% Select those bars outside the cone boundary (both ends bars are outside of the cone boundary) 
TH2BARS                    = reshape(GCOORD_SPH_c(1,bar2node_els_in_cone_iso),[],2);
theta_bars_out_cone        = sum(TH2BARS > theta_cone*pi/180 & TH2BARS < (180-theta_cone)*pi/180,2); % boolean vector for bars having both ends outside the cone boundary
                                    % 0 --> both ends bar are inside the cone boundary 
                                    % 1 --> one end bar is inside the cone boundary and the other one is outside the cone boundary 
                                    % 2 --> both ends bar are outside the cone boundary 
bar2node_els_in_cone_iso_with_2_bars_out_cone = bar2node_els_in_cone_iso(theta_bars_out_cone == 2,:);

% Remove these ouside boundary cone bars in the bar2node_els_in_cone list 
[~,LOCA]  = ismember(bar2node_els_in_cone_iso_with_2_bars_out_cone,bar2node_els_in_cone,'rows'); % compute positions of bar2node_els_in_cone_iso_with_2_bars_out_cone
                                                                                                 % in the bar2node_els_in_cone list
[~,LOCA2] = ismember(fliplr(bar2node_els_in_cone_iso_with_2_bars_out_cone),bar2node_els_in_cone,'rows'); % this is just in case some bars were flipped 
                                                                                                         % when using find_unique_bars subfunction
LOCA(LOCA == 0)              = LOCA2(LOCA2 ~= 0); clear LOCA2;
bar2node_els_in_cone(LOCA,:) = []; % remove these ouside boundary cone bars in the bar2node_els_in_cone list 

% Coordinates of the nodes defining each bar's end points
thetaBarEnds_els_in_cone   = reshape(GCOORD_SPH_c_rot(1,bar2node_els_in_cone'),2,[]);
phiBarEnds_els_in_cone     = reshape(GCOORD_SPH_c_rot(2,bar2node_els_in_cone'),2,[]);
rBarEnds_els_in_cone       = reshape(GCOORD_SPH_c_rot(3,bar2node_els_in_cone'),2,[]);

% Create new node at each bar mid-point
thetaBarMids_els_in_cone   = 0.5*sum(thetaBarEnds_els_in_cone,1);
phiBarMids_els_in_cone     = 0.5*sum(phiBarEnds_els_in_cone,1);
rBarMids_els_in_cone       = 0.5*sum(rBarEnds_els_in_cone,1);
GCOORD_SPH_els_in_cone_rot = [thetaBarMids_els_in_cone; ...
                              phiBarMids_els_in_cone;   ...
                              rBarMids_els_in_cone];
 
% Convert to Cartesian coordinates
GCOORD_els_in_cone_rot     = spherical2cartesian(GCOORD_SPH_els_in_cone_rot);

% Undo the rotation (rotate 90° around X axis clockwise)
RR_X_90_CW                 = [ 1  0  0 ; ...
                               0  0  1 ; ...
                               0 -1  0 ];
GCOORD_els_in_cone         = RR_X_90_CW * GCOORD_els_in_cone_rot;

% Convert to spherical coodinates
GCOORD_SPH_els_in_cone     = cartesian2spherical(GCOORD_els_in_cone);

end % END OF SUBFUNCTION mid_side_nodes_sph_rot_frame_els_in_cone

% #####################################################################################################

function [GCOORD_curved,GCOORD_SPH,EL2NOD] = mid_side_nodes_sph_frame_all_els...
    (GCOORD_SPH_c,EL2NOD_c,GCOORD_SPH_els_in_cone,bar2node_els_in_cone,GCOORD_SPH_els_out_cone_cross_2pi,bar2node_els_out_cone_cross_2pi)
% Usage: [GCOORD_curved,GCOORD_SPH,EL2NOD] = mid_side_nodes_sph_frame_all_els...
%   (GCOORD_SPH_c,EL2NOD_c,GCOORD_SPH_els_in_cone,bar2node_els_in_cone,GCOORD_SPH_els_out_cone_cross_2pi,bar2node_els_cross_2pi_out_cone)
%
% Purpose: 
%   Compute mid-side nodes for all elements in the original sph frame using
%   mid-side nodes computed in 'mid_side_nodes_sph_rot_frame_els_in_cone'
%   and 'mid_side_nodes_sph_rot_frame_els_cross_2pi_out_cone' subfunctions.
%   The steps are:  
%       - Compute bars for all elements (bar2node).
%       - Compute mid-side nodes in straight edges of the original
%         spherical frame.
%       - Find bar positions of bar2node_els_in_cone in bar2node.
%       - Substitute GCOORD_SPH_els_in_cone in the mid-side nodes computed
%         in straight edges of the original spherical frame using the bar
%         positions of bar2node_els_in_cone in bar2node.
%       - Find bar positions of bar2node_els_out_cone_cross_2pi in bar2node
%       - Substitute GCOORD_SPH_els_out_cone_cross_2pi in the mid-side
%         nodes computed in straight edges of the original spherical frame
%         using the bar positions of bar2node_els_out_cone_cross_2pi in
%         bar2node.
%
% Input:
%   GCOORD_SPH_c                      : [matrix] : spherical coordinates of
%                                                  all nodes (4 nodel)
%                                                  (3 x nnod_c)
%   EL2NOD_c                          : [matrix] : connectivity matrix 
%                                                  (4 x nel)
%   GCOORD_SPH_els_in_cone            : [matrix] : spherical coordinates of
%                                                  nodes belonging those
%                                                  elements in the cone
%                                                  (els_in_cone = 
%                                                  els_in_cone_no_iso +
%                                                  els_in_cone_iso)
%   bar2node_els_in_cone              : [matrix] : bar list for elements in
%                                                  the cone (els_in_cone =
%                                                  els_in_cone_no_iso +
%                                                  els_in_cone_iso)
%   GCOORD_SPH_els_out_cone_cross_2pi : [matrix] : spherical coordinates of
%                                                  nodes belonging those
%                                                  elements outside the
%                                                  cone and crossing 
%                                                  phi = 2pi
%   bar2node_els_out_cone_cross_2pi   : [matrix] : bar list for elements
%                                                  outside the cone and
%                                                  crossing phi = 2pi
%
% Output:
%   GCOORD_curved : [matrix] : Cartesian coordinates of all nodes new mesh
%                              (3 x ?) 
%   GCOORD_SPH    : [matrix] : spherical coordinates of all nodes new mesh
%                              (3 x ?) 
%   EL2NOD        : [matrix] : connectivity matrix (10 x nel)
%
% JMT Aug 2016

% Size of coarse mesh
nel_c              = size(EL2NOD_c,2);     % number of elements
nnod_c             = size(GCOORD_SPH_c,2); % number of nodes

% Create a pointer bars defined by their end-nodes (e.g. bar 1 has end-nodes [1 5], bar 2 has [5 2],...)
bar2node           = reshape( EL2NOD_c([ 1  2 ... % 1st edge of parent
                                         2  3 ... % 2nd edge of parent
                                         3  4 ... % 3rd edge of parent
                                         4  1 ... % 4th edge of parent
                                         1  3 ... % 5th edge of parent
                                         4  2 ... % 6th edge of parent
                                             ]',:),2,[] )';

% Find the bars that are shared by neighboring elements and return a unique list of bars
[bar2node,~,ib]    = find_unique_bars(bar2node); % *SUBFUNCTION*
nbars              = size(bar2node,1);    % number of unique bars
EL2BAR             = reshape(ib,6,nel_c); % element to bar connectivity after doubles have been merged

% Coordinates of the nodes defining each bar's end points
thetaBarEnds       = reshape(GCOORD_SPH_c(1,bar2node'),2,[]);
phiBarEnds         = reshape(GCOORD_SPH_c(2,bar2node'),2,[]);
rBarEnds           = reshape(GCOORD_SPH_c(3,bar2node'),2,[]);

% Check if phi angles (for bar nodes) are separated an angular distance bigger than pi 
phiBarEnds         = check_ang_dist_phi(GCOORD_SPH_c,bar2node,phiBarEnds);

% Create new node at each bar mid-point
thetaBarMids       = 0.5*sum(thetaBarEnds,1);
phiBarMids         = 0.5*sum(phiBarEnds,1);
rBarMids           = 0.5*sum(rBarEnds,1);

% Compute positions of bar2node_els_in_cone in the bar2node list
[~,LOCB]           = ismember(bar2node_els_in_cone,bar2node,'rows');
[~,LOCB2]          = ismember(fliplr(bar2node_els_in_cone),bar2node,'rows'); % this is just in case some bars were flipped when using find_unique_bars subfunction
LOCB(LOCB == 0)    = LOCB2(LOCB2 ~= 0); clear LOCB2;
% Substitute mid-side nodes of elements in the cone (elements within the cone boundary + elements crossing the cone boundary) 
thetaBarMids(LOCB) = GCOORD_SPH_els_in_cone(1,:);
phiBarMids(LOCB)   = GCOORD_SPH_els_in_cone(2,:);
rBarMids(LOCB)     = GCOORD_SPH_els_in_cone(3,:);

% Compute positions of bar2node_els_out_cone_cross_2pi in the bar2node list
[~,LOCD]           = ismember(bar2node_els_out_cone_cross_2pi,bar2node,'rows');
[~,LOCD2]          = ismember(fliplr(bar2node_els_out_cone_cross_2pi),bar2node,'rows'); % this is just in case some bars were flipped when using find_unique_bars subfunction
LOCD(LOCD == 0)    = LOCD2(LOCD2 ~= 0); clear LOCD2;
% Substitute mid-side nodes of elements outside the cone and crossing phi = 2pi 
thetaBarMids(LOCD) = GCOORD_SPH_els_out_cone_cross_2pi(1,:);
phiBarMids(LOCD)   = GCOORD_SPH_els_out_cone_cross_2pi(2,:);
rBarMids(LOCD)     = GCOORD_SPH_els_out_cone_cross_2pi(3,:);

% Storage for quadratic order (10-node) connectivity matrix
EL2NOD             = zeros(10,nel_c,'uint32'); 
EL2NOD(1: 4,:)     = EL2NOD_c;
EL2NOD(5:10,:)     = nnod_c+EL2BAR;

% Storage for quadratic order (10-node) node coordinates
nnod                        = nnod_c + nbars;
GCOORD_SPH                  = zeros(3,nnod);
GCOORD_SPH(:,1:nnod_c)      = GCOORD_SPH_c;
GCOORD_SPH(1,nnod_c+1:nnod) = thetaBarMids;
GCOORD_SPH(2,nnod_c+1:nnod) = phiBarMids;
GCOORD_SPH(3,nnod_c+1:nnod) = rBarMids;

% Convert to Cartesian coordinates
GCOORD_curved               = spherical2cartesian(GCOORD_SPH);

end % END OF SUBFUNCTION mid_side_nodes_sph_frame_all_els

% #####################################################################################################

function distortion = plot_distortion_mid_side_nodes...
    (GCOORD_c,GCOORD_SPH_c,EL2NOD_c,EL2NOD,els_in_cone_iso)
% Usage: distortion = plot_distortion_mid_side_nodes...
%   (GCOORD_c,GCOORD_SPH_c,EL2NOD_c,EL2NOD,els_in_cone_iso)
%
% Purpose: 
%   Compute mid-side distortion (only for debugging purposes) for those
%   elements crossing the cone boundary and having at least 1 edge outside
%   the cone (elements with curved spherical edges in the rotated spherical
%   frame).
%   Firstly we compute the distance between:
%       - mid-side nodes computed in the original spherical coordinate
%         system (GCOORD_SPH) 
%       - mid-side nodes computed in the rotated spherical coordinate
%         system (GCOORD_SPH_rot) 
%   In order to subtract both positions, we need to do a 90° rotation of
%   the Cartesian coordinates of mid-side nodes computed in the original
%   spherical coordinate system, and then, transform to spherical
%   coordinates (GCOORD_SPH_rot_to_compare).
%   Then we divide this distance by the length of each edge where those
%   mid-side nodes belong. The length of each edge is computed as the sum
%   of distances from a vertex to mid-side + mid-side to other vertex.
%   Finelly we multiply by 100 to give a distortion percentage
%
% Input:
%   GCOORD_c        : [matrix] : Cartesian coordinates of all nodes 
%                                (4 nodel) (3 x nnod_c)
%   GCOORD_SPH_c    : [matrix] : spherical coordinates of all nodes 
%                                (4 nodel) (3 x nnod_c)
%   EL2NOD_c        : [matrix] : connectivity matrix (4 x nel)
%   EL2NOD          : [matrix] : connectivity matrix (10 x nel)
%   els_in_cone_iso : [vector] : indices for those elements crossing the
%                                cone boundary and having at least 1 edge
%                                (bar) outside the cone boundary
%
% Output:
%   distortion      : [matrix] : distortion in % for each mid-side node of
%                                each element crossing the cone boundary
%                                and having at least 1 edge outside the
%                                cone (elements with curved spherical edges
%                                in the rotated spherical frame). 
%
% JMT Aug 2016

% ==============================================================================================================
% COMPUTE MIDPOINTS IN SPHERICAL COORDINATES (CURVED EDGES) FOR ALL ELEMENTS  
% ==============================================================================================================
% Size of coarse mesh
nnod_c                            = size(GCOORD_SPH_c,2); % number of nodes

% Create a pointer bars defined by their end-nodes (e.g. bar 1 has end-nodes [1 5], bar 2 has [5 2],...)
bar2node                          = reshape( EL2NOD_c([ 1  2 ... % 1st edge of parent
                                                        2  3 ... % 2nd edge of parent
                                                        3  4 ... % 3rd edge of parent
                                                        4  1 ... % 4th edge of parent
                                                        1  3 ... % 5th edge of parent
                                                        4  2 ... % 6th edge of parent
                                                            ]',:),2,[] )';

% Find the bars that are shared by neighboring elements and return a unique list of bars
[bar2node,~,~]                    = find_unique_bars(bar2node); % *SUBFUNCTION*
nbars                             = size(bar2node,1); % number of unique bars

% Coordinates of the nodes defining each bar's end points
thetaBarEnds                      = reshape(GCOORD_SPH_c(1,bar2node'),2,[]);
phiBarEnds                        = reshape(GCOORD_SPH_c(2,bar2node'),2,[]);
rBarEnds                          = reshape(GCOORD_SPH_c(3,bar2node'),2,[]);

% Check if phi angles (for bar nodes) are separated an angular distance bigger than pi 
phiBarEnds                        = check_ang_dist_phi(GCOORD_SPH_c,bar2node,phiBarEnds);

% Create new node at each bar mid-point
thetaBarMids_no_rot               = 0.5*sum(thetaBarEnds,1);
phiBarMids_no_rot                 = 0.5*sum(phiBarEnds,1);
rBarMids_no_rot                   = 0.5*sum(rBarEnds,1);

% Check if any new node in phiBarMids has phi > 2pi, if so subtract 2pi
phiBarMids_no_rot(phiBarMids_no_rot > 2*pi) = phiBarMids_no_rot(phiBarMids_no_rot > 2*pi) - 2*pi;

% Storage for quadratic order (10-node) node coordinates
nnod                               = nnod_c + nbars;
GCOORD_SPH_no_rot                  = zeros(3,nnod);
GCOORD_SPH_no_rot(:,1:nnod_c)      = GCOORD_SPH_c;
GCOORD_SPH_no_rot(1,nnod_c+1:nnod) = thetaBarMids_no_rot;
GCOORD_SPH_no_rot(2,nnod_c+1:nnod) = phiBarMids_no_rot;
GCOORD_SPH_no_rot(3,nnod_c+1:nnod) = rBarMids_no_rot;

GCOORD_curved_no_rot               = spherical2cartesian(GCOORD_SPH_no_rot);

% ==============================================================================================================
% COMPUTE MIDPOINTS IN SPHERICAL COORDINATES (CURVED EDGES) FOR ALL ELEMENTS AFTER A 90° ROTATION AROUND X AXIS
% ==============================================================================================================
% Make a 90° counterclockwise rotation around the X axis of the elements in the cone
RR_X_90_CCW                        = [ 1  0  0 ; ...
                                       0  0 -1 ; ...
                                       0  1  0 ];
GCOORD_c_rot                       = RR_X_90_CCW * GCOORD_c;

% Convert to spherical coordinates
GCOORD_SPH_c_rot                   = cartesian2spherical(GCOORD_c_rot);

% Coordinates of the nodes defining each bar's end points
thetaBarEnds_rot                   = reshape(GCOORD_SPH_c_rot(1,bar2node'),2,[]);
phiBarEnds_rot                     = reshape(GCOORD_SPH_c_rot(2,bar2node'),2,[]);
rBarEnds_rot                       = reshape(GCOORD_SPH_c_rot(3,bar2node'),2,[]);

% Check if phi angles (for bar nodes) are separated an angular distance bigger than pi
phiBarEnds_rot                     = check_ang_dist_phi(GCOORD_SPH_c_rot,bar2node,phiBarEnds_rot);

% Create new node at each bar mid-point
thetaBarMids_rot                   = 0.5*sum(thetaBarEnds_rot,1);
phiBarMids_rot                     = 0.5*sum(phiBarEnds_rot,1);
rBarMids_rot                       = 0.5*sum(rBarEnds_rot,1);

% Check if any new node in phiBarMids has phi > 2pi, if so subtract 2pi
phiBarMids_rot(phiBarMids_rot > 2*pi) = phiBarMids_rot(phiBarMids_rot > 2*pi) - 2*pi;

% Storage for quadratic order (10-node) node coordinates
GCOORD_SPH_rot                     = zeros(3,nnod);
GCOORD_SPH_rot(:,1:nnod_c)         = GCOORD_SPH_c_rot;
GCOORD_SPH_rot(1,nnod_c+1:nnod)    = thetaBarMids_rot;
GCOORD_SPH_rot(2,nnod_c+1:nnod)    = phiBarMids_rot;
GCOORD_SPH_rot(3,nnod_c+1:nnod)    = rBarMids_rot;

GCOORD_curved_rot                  = spherical2cartesian(GCOORD_SPH_rot);

% Undo the rotation
RR_X_90_CW                         = [ 1  0  0 ; ...
                                       0  0  1 ; ...
                                       0 -1  0 ];
GCOORD_curved_undo_rot             = RR_X_90_CW * GCOORD_curved_rot;

% ==============================================================================================================
% COMPUTE DISTORTION
% ==============================================================================================================
distortion = zeros(size(els_in_cone_iso,2),6);
nodV       = [1  2;
              2  3;
              3  4;
              4  1;
              1  3;
              4  2];
for i = 1:size(els_in_cone_iso,2)
    
    mid_side_nodes_original_sph_frame = GCOORD_curved_no_rot(:,EL2NOD(5:10,els_in_cone_iso(i)));
    mid_side_nodes_rotated_sph_frame  = GCOORD_curved_undo_rot(:,EL2NOD(5:10,els_in_cone_iso(i)));
    
    d_mid_side_nodes = sqrt((mid_side_nodes_original_sph_frame(1,:)-mid_side_nodes_rotated_sph_frame(1,:)).^2 + ...
                            (mid_side_nodes_original_sph_frame(2,:)-mid_side_nodes_rotated_sph_frame(2,:)).^2 + ...
                            (mid_side_nodes_original_sph_frame(3,:)-mid_side_nodes_rotated_sph_frame(3,:)).^2);
    
    vertices         = GCOORD_curved_undo_rot(:,EL2NOD(1:4,els_in_cone_iso(i)));
    d_vertices1      = sqrt((vertices(1,nodV(:,1))-mid_side_nodes_rotated_sph_frame(1,:)).^2 + ...
                            (vertices(2,nodV(:,1))-mid_side_nodes_rotated_sph_frame(2,:)).^2 + ...
                            (vertices(3,nodV(:,1))-mid_side_nodes_rotated_sph_frame(3,:)).^2);
    d_vertices2      = sqrt((vertices(1,nodV(:,2))-mid_side_nodes_rotated_sph_frame(1,:)).^2 + ...
                            (vertices(2,nodV(:,2))-mid_side_nodes_rotated_sph_frame(2,:)).^2 + ...
                            (vertices(3,nodV(:,2))-mid_side_nodes_rotated_sph_frame(3,:)).^2);
    d_vertices       = d_vertices1 + d_vertices2;
    distortion(i,:)  = (d_mid_side_nodes./d_vertices)*100;
    
end

fprintf('Max distortion in mid-side nodes = %4.2f %%\n',max(max(distortion)));

% ==============================================================================================================
% PLOT
% ==============================================================================================================
figure(81)
clf
x      = 1:size(els_in_cone_iso,2);
y      = max(distortion,[],2);
y_mean = mean(y)*ones(size(els_in_cone_iso,2),1);
plot(x,y,'+-')
hold on
plot(x,y_mean,'r-')
legend('distortion','distortion mean','Location','SouthEast')
xlabel('Elements')
ylabel('Distortion')
title('Maximun distortion for each element with curved edges in the rotated spherical frame')

end % END OF SUBFUNCTION plot_distortion_mid_side_nodes

% #####################################################################################################

function plot_figures(GCOORD_c,GCOORD_curved,GCOORD_SPH,EL2NOD,els_in_cone_iso,els_within_cone,els_cone_bnd,els_out_cone,els_out_cone_cross_2pi,distortion,theta_cone)

[el,~] = find(distortion == max(max(distortion)));

el_in_cone_iso_plot     = els_in_cone_iso(el(1)); % index of element on the cone boundary to be plotted

[iI,~]                  = ismember(EL2NOD(1:4,els_within_cone),EL2NOD(1:4,el_in_cone_iso_plot));
i_els_within_cone_attached_el_bnd_cone  = find(sum(iI)==3); % indices for elements within the cone sharing 3 vertices with element on the cone bnd
if isempty(i_els_within_cone_attached_el_bnd_cone)
    i_els_within_cone_attached_el_bnd_cone = find(sum(iI)==2); % indices for elements within the cone sharing 2 vertices with element on the cone bnd
    if isempty(i_els_within_cone_attached_el_bnd_cone)
        i_els_within_cone_attached_el_bnd_cone = find(sum(iI)==1); % indices for elements within the cone sharing 1 vertex with element on the cone bnd
    end
end
el_within_cone_plot     = els_within_cone(i_els_within_cone_attached_el_bnd_cone(1)); % index of element within the cone to be plotted

[iO,~]                  = ismember(EL2NOD(1:4,els_out_cone),EL2NOD(1:4,el_in_cone_iso_plot));
i_els_outside_cone_attached_el_bnd_con = find(sum(iO)==3); % indices for element outside the cone sharing 3 vertices with element on the cone bnd
if isempty(i_els_outside_cone_attached_el_bnd_con)
    i_els_outside_cone_attached_el_bnd_con = find(sum(iO)==2); % indices for element outside the cone sharing 2 vertices with element on the cone bnd
    if isempty(i_els_outside_cone_attached_el_bnd_con)
        i_els_outside_cone_attached_el_bnd_con = find(sum(iO)==1); % indices for element outside the cone sharing 1 vertex with element on the cone bnd
    end
end
el_out_cone_plot           = els_out_cone(i_els_outside_cone_attached_el_bnd_con(1)); % index of element outside the cone to be plotted
el_out_cone_cross_2pi_plot = els_out_cone_cross_2pi(1);


figure(82)
clf
GCOORD_els_cone_bnd_vertices = unique(GCOORD_SPH(:,EL2NOD(1:4,[el_in_cone_iso_plot el_within_cone_plot el_out_cone_plot]))','rows','stable');
GCOORD_els_cone_bnd_vertices = GCOORD_els_cone_bnd_vertices';
GCOORD_els_cone_bnd_vertices(3,:) = GCOORD_els_cone_bnd_vertices(3,:)/1000;
GCOORD_els_cone_bnd          = unique(GCOORD_SPH(:,EL2NOD(:,[el_in_cone_iso_plot el_within_cone_plot el_out_cone_plot]))','rows','stable');
GCOORD_els_cone_bnd          = GCOORD_els_cone_bnd';
GCOORD_els_cone_bnd(3,:)     = GCOORD_els_cone_bnd(3,:)/1000;
GCOORD_SPH_plot              = GCOORD_SPH;
GCOORD_SPH_plot(3,:)         = GCOORD_SPH_plot(3,:)/1000;
lightMag  = 0.80*[1 0.7 1]; % colour for elements within the cone
lightYell = 0.80*[1 1 0];   % colour for elements crossing the cone boundary
lightCyan = 0.80*[0 1 1];   % colour for elements crossing the cone boundary and having at least 1 edge outside the cone boundary 
                            % (elements with curved spherical sides in the spherical rotated frame)
tetramesh(EL2NOD(:,el_within_cone_plot)',GCOORD_SPH_plot','FaceColor',lightMag,'FaceAlpha',0.3)
hold on
tetramesh(EL2NOD(:,el_in_cone_iso_plot)',GCOORD_SPH_plot','FaceColor',lightCyan,'FaceAlpha',0.3)
tetramesh(EL2NOD(:,el_out_cone_plot)',GCOORD_SPH_plot','FaceColor',[1 1 1],'FaceAlpha',0)
scatter3(GCOORD_els_cone_bnd(1,:),GCOORD_els_cone_bnd(2,:),GCOORD_els_cone_bnd(3,:), ...
    'MarkerEdgeColor','k','MarkerFaceColor',[0 1 0])
scatter3(GCOORD_els_cone_bnd_vertices(1,:),GCOORD_els_cone_bnd_vertices(2,:),GCOORD_els_cone_bnd_vertices(3,:), ...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
% axis([0 pi 0 2*pi 3 6.5])
view(142.5,30)
grid on
xlabel('\theta (rad)')
ylabel('\phi (rad)')
zlabel('r (10^3 km)')
title('Spherical coordinates')
legend('Element within the cone bnd (Straight edges in the rotated spherical frame)', ...
       'Element crossing the cone bnd (Curved edges in the rotated spherical frame)', ...
       'Element outside the cone bnd (Straight edges in the original spherical frame)', ...
       'Mid-side node','Vertices','Location','SouthOutside')

figure(83)
clf
tetramesh(EL2NOD(:,el_within_cone_plot)',GCOORD_curved','FaceColor',lightMag,'FaceAlpha',0.3)
hold on
tetramesh(EL2NOD(:,el_in_cone_iso_plot)',GCOORD_curved','FaceColor',lightCyan,'FaceAlpha',0.3)
tetramesh(EL2NOD(:,el_out_cone_plot)',GCOORD_curved','FaceColor',[1 1 1],'FaceAlpha',0)
GCOORD_els_cone_bnd_vertices = unique(GCOORD_c(:,EL2NOD(1:4,[el_in_cone_iso_plot el_within_cone_plot el_out_cone_plot]))','rows','stable');
GCOORD_els_cone_bnd_vertices = GCOORD_els_cone_bnd_vertices';
GCOORD_els_cone_bnd          = unique(GCOORD_curved(:,EL2NOD(:,[el_in_cone_iso_plot el_within_cone_plot el_out_cone_plot]))','rows','stable');
GCOORD_els_cone_bnd          = GCOORD_els_cone_bnd';
scatter3(GCOORD_els_cone_bnd(1,:),GCOORD_els_cone_bnd(2,:),GCOORD_els_cone_bnd(3,:), ...
    'MarkerEdgeColor','k','MarkerFaceColor',[0 1 0])
scatter3(GCOORD_els_cone_bnd_vertices(1,:),GCOORD_els_cone_bnd_vertices(2,:),GCOORD_els_cone_bnd_vertices(3,:), ...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
view(142.5,30)
grid on
title('Cartesian coordinates')
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
legend('Element within the cone bnd (Straight edges in the rotated spherical frame)', ...
       'Element crossing the cone bnd (Curved edges in the rotated spherical frame)', ...
       'Element outside the cone bnd (Straight edges in the original spherical frame)', ...
       'Mid-side node','Vertices','Location','SouthOutside')

figure(84)
clf
th    = theta_cone*pi/180;
h     = 6371;
r     = h*tan(th);
[R,A] = meshgrid(linspace(0,r,2),linspace(0,2*pi,21));
X     = R .* cos(A);
Y     = R .* sin(A);
Z     = R/tan(th);
% Cone around the z-axis, point at the origin
surface(X,Y,Z,'FaceColor','none','EdgeColor',lightMag)
hold on
surface(X,Y,-Z,'FaceColor','none','EdgeColor',lightMag)
lightGrey                    = 0.90*[1 1 1];
[x_sphere,y_sphere,z_sphere] = sphere(20);
x_sphere                     = x_sphere*3471;
y_sphere                     = y_sphere*3471;
z_sphere                     = z_sphere*3471;
axis equal
surface(x_sphere,y_sphere,z_sphere,'FaceColor','none','EdgeColor',lightGrey)
hold on
[x_sphere,y_sphere,z_sphere] = sphere(30);
x_sphere                     = x_sphere*6371;
y_sphere                     = y_sphere*6371;
z_sphere                     = z_sphere*6371;
surface(x_sphere,y_sphere,z_sphere,'FaceColor','none','EdgeColor',lightGrey)
axis([-6371 6371 -6371 6371 -6371 6371])
view(142.5,30)
grid on
tetramesh(EL2NOD(:,el_within_cone_plot)',GCOORD_curved','FaceColor',lightMag,'FaceAlpha',0.3)
hold on
tetramesh(EL2NOD(:,el_in_cone_iso_plot)',GCOORD_curved','FaceColor',lightCyan,'FaceAlpha',0.3)
tetramesh(EL2NOD(:,el_out_cone_plot)',GCOORD_curved','FaceColor',[1 1 1],'FaceAlpha',0)
tetramesh(EL2NOD(:,el_out_cone_cross_2pi_plot)',GCOORD_curved','FaceColor',lightYell,'FaceAlpha',0.3)
GCOORD_els_cone_bnd_vertices = unique(GCOORD_c(:,EL2NOD(1:4,[el_in_cone_iso_plot el_within_cone_plot el_out_cone_plot el_out_cone_cross_2pi_plot]))','rows','stable');
GCOORD_els_cone_bnd_vertices = GCOORD_els_cone_bnd_vertices';
GCOORD_els_cone_bnd          = unique(GCOORD_curved(:,EL2NOD(:,[el_in_cone_iso_plot el_within_cone_plot el_out_cone_plot el_out_cone_cross_2pi_plot]))','rows','stable');
GCOORD_els_cone_bnd          = GCOORD_els_cone_bnd';
scatter3(GCOORD_els_cone_bnd(1,:),GCOORD_els_cone_bnd(2,:),GCOORD_els_cone_bnd(3,:), ...
    'MarkerEdgeColor','k','MarkerFaceColor',[0 1 0])
scatter3(GCOORD_els_cone_bnd_vertices(1,:),GCOORD_els_cone_bnd_vertices(2,:),GCOORD_els_cone_bnd_vertices(3,:), ...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
title('Elements location within the spherical shell')
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')

figure(85)
clf
subplot(1,3,1)
% Cone around the z-axis, point at the origin
surface(X,Y,Z,'FaceColor','none','EdgeColor',lightMag)
hold on
surface(X,Y,-Z,'FaceColor','none','EdgeColor',lightMag)
[x_sphere,y_sphere,z_sphere] = sphere(20);
x_sphere = x_sphere*3471;
y_sphere = y_sphere*3471;
z_sphere = z_sphere*3471;
axis equal
surface(x_sphere,y_sphere,z_sphere,'FaceColor','none','EdgeColor',lightGrey)
hold on
[x_sphere,y_sphere,z_sphere] = sphere(30);
x_sphere = x_sphere*6371;
y_sphere = y_sphere*6371;
z_sphere = z_sphere*6371;
axis equal
surface(x_sphere,y_sphere,z_sphere,'FaceColor','none','EdgeColor',lightGrey)
hold on
axis([-6371 6371 -6371 6371 -6371 6371])
view(142.5,30)
grid on
tetramesh(EL2NOD(:,els_within_cone)',GCOORD_curved','FaceColor',lightMag,'FaceAlpha',0.3)
GCOORD_els_inside_cone_vertices = unique(GCOORD_c(:,EL2NOD(1:4,els_within_cone))','rows','stable');
GCOORD_els_inside_cone_vertices = GCOORD_els_inside_cone_vertices';
GCOORD_els_in_cone              = unique(GCOORD_curved(:,EL2NOD(:,els_within_cone))','rows','stable');
GCOORD_els_in_cone              = GCOORD_els_in_cone';
scatter3(GCOORD_els_in_cone(1,:),GCOORD_els_in_cone(2,:),GCOORD_els_in_cone(3,:), ...
    'MarkerEdgeColor','k','MarkerFaceColor',[0 1 0])
scatter3(GCOORD_els_inside_cone_vertices(1,:),GCOORD_els_inside_cone_vertices(2,:),GCOORD_els_inside_cone_vertices(3,:), ...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
title('Elements within the cone (Straight edges in the rotated spherical frame)')

subplot(1,3,2)
% Cone around the z-axis, point at the origin
surface(X,Y,Z,'FaceColor','none','EdgeColor',lightMag)
hold on
surface(X,Y,-Z,'FaceColor','none','EdgeColor',lightMag)
[x_sphere,y_sphere,z_sphere] = sphere(20);
x_sphere = x_sphere*3471;
y_sphere = y_sphere*3471;
z_sphere = z_sphere*3471;
axis equal
surface(x_sphere,y_sphere,z_sphere,'FaceColor','none','EdgeColor',lightGrey)
hold on
[x_sphere,y_sphere,z_sphere] = sphere(30);
x_sphere = x_sphere*6371;
y_sphere = y_sphere*6371;
z_sphere = z_sphere*6371;
axis equal
surface(x_sphere,y_sphere,z_sphere,'FaceColor','none','EdgeColor',lightGrey)
hold on
axis([-6371 6371 -6371 6371 -6371 6371])
view(142.5,30)
grid on
tetramesh(EL2NOD(:,els_cone_bnd)',GCOORD_curved','FaceColor',lightYell,'FaceAlpha',0.3)
GCOORD_els_inside_cone_vertices = unique(GCOORD_c(:,EL2NOD(1:4,els_cone_bnd))','rows','stable');
GCOORD_els_inside_cone_vertices = GCOORD_els_inside_cone_vertices';
GCOORD_els_in_cone              = unique(GCOORD_curved(:,EL2NOD(:,els_cone_bnd))','rows','stable');
GCOORD_els_in_cone              = GCOORD_els_in_cone';
scatter3(GCOORD_els_in_cone(1,:),GCOORD_els_in_cone(2,:),GCOORD_els_in_cone(3,:), ...
    'MarkerEdgeColor','k','MarkerFaceColor',[0 1 0])
scatter3(GCOORD_els_inside_cone_vertices(1,:),GCOORD_els_inside_cone_vertices(2,:),GCOORD_els_inside_cone_vertices(3,:), ...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
title('Elements crossing the cone boundary')

subplot(1,3,3)
% Cone around the z-axis, point at the origin
surface(X,Y,Z,'FaceColor','none','EdgeColor',lightMag)
hold on
surface(X,Y,-Z,'FaceColor','none','EdgeColor',lightMag)
[x_sphere,y_sphere,z_sphere] = sphere(20);
x_sphere = x_sphere*3471;
y_sphere = y_sphere*3471;
z_sphere = z_sphere*3471;
axis equal
surface(x_sphere,y_sphere,z_sphere,'FaceColor','none','EdgeColor',lightGrey)
hold on
[x_sphere,y_sphere,z_sphere] = sphere(30);
x_sphere = x_sphere*6371;
y_sphere = y_sphere*6371;
z_sphere = z_sphere*6371;
axis equal
surface(x_sphere,y_sphere,z_sphere,'FaceColor','none','EdgeColor',lightGrey)
hold on
axis([-6371 6371 -6371 6371 -6371 6371])
view(142.5,30)
grid on
tetramesh(EL2NOD(:,els_in_cone_iso)',GCOORD_curved','FaceColor',lightCyan,'FaceAlpha',0.3)
GCOORD_els_inside_cone_vertices = unique(GCOORD_c(:,EL2NOD(1:4,els_in_cone_iso))','rows','stable');
GCOORD_els_inside_cone_vertices = GCOORD_els_inside_cone_vertices';
GCOORD_els_in_cone              = unique(GCOORD_curved(:,EL2NOD(:,els_in_cone_iso))','rows','stable');
GCOORD_els_in_cone              = GCOORD_els_in_cone';
scatter3(GCOORD_els_in_cone(1,:),GCOORD_els_in_cone(2,:),GCOORD_els_in_cone(3,:), ...
    'MarkerEdgeColor','k','MarkerFaceColor',[0 1 0])
scatter3(GCOORD_els_inside_cone_vertices(1,:),GCOORD_els_inside_cone_vertices(2,:),GCOORD_els_inside_cone_vertices(3,:), ...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
title('Elements crossing the cone boundary having at least 1 edge outside the cone boundary (Curved edges in the rotated spherical frame)')

end % END OF SUBFUNCTION plot_figures

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