function [els,gTH_PT,OPTS_search] = tsearch2_sph(GCOORD_SPH,EL2NOD,gTH_PT,OPTS_search,MESH)
% Usage: [els,gTH_PT,OPTS_search] = tsearch2_sph(GCOORD_SPH,EL2NOD,gTH_PT,OPTS_search,MESH)
%
% Purpose: 
%   Locate points with coordinates "gTH_PT" in a 3D FE spherical mesh 
%
% Input:
%   GCOORD_SPH          : [matrix]    : spherical coodinates (theta,phi,r) 
%                                       (in rad ,rad and km) of the mesh
%   EL2NOD              : [matrix]    : connectivity matrix for GCOORD_SPH
%   gTH_PT              : [matrix]    : coordinates of points to be located
%   OPTS_search         : [structure] : options for tsearch2
%   MESH                : [structure] : FE mesh parameters
%
% Output:
%   els                 : [rowvector] : element in which each point is 
%                                       located 
%   gTH_PT              : [matrix]    : coordinates of points to be located
%   OPTS_search         : [structure] : options for tsearch2_sph
%
% Part of M3TRI 2D convection code family,
% developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de
% For numerical methods see online Ph.D. thesis
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)

% JH Feb 2011
% JH Apr 2015 : added tsearch2
% JMT May 2016 : Now it works for polar-coordinate meshes
% JMT Oct 2016 : works for spherical coordinates

%===========================================================================================
% LOAD WS AND COMPUTE INFO RELATED TO THE CONE -> SPLIT POINTS TO BE LOCATED IN TWO GROUPS:
%   - POINTS TO BE LOCATED THAT ARE INSIDE THE CONE
%   - POINTS TO BE LOCATED THAT ARE OUTSIDE THE CONE
%===========================================================================================
els_out_cone_no_cross_2pi  = MESH.els_out_cone_no_cross_2pi;
els_out_cone_cross_2pi     = MESH.els_out_cone_cross_2pi;
els_in_cone_no_iso         = MESH.els_in_cone_no_iso;
els_in_cone_iso            = MESH.els_in_cone_iso;
WS_els_in_cone_rot_X_90    = OPTS_search.WS_els_in_cone_rot_X_90;
WS_els_out_cone            = OPTS_search.WS_els_out_cone;
WS_els_cross_2pi_rot_Z_180 = OPTS_search.WS_els_cross_2pi_rot_Z_180;
theta_cone                 = MESH.theta_cone*pi/180; % theta cone in rad
theta_north                = theta_cone;             % theta angle from the +Z axis to the generatrix
theta_south                = pi-theta_cone;          % theta angle from the -Z axis to the generatrix
points_in_cone             = gTH_PT(1,:) <= theta_north | gTH_PT(1,:) >= theta_south; % boolean vector for nodes in the cone
points_out_cone            = find(~points_in_cone);
nnod                       = max(max(EL2NOD(1:4,:)));
els                        = zeros(size(gTH_PT,2),4);
GCOORD                     = spherical2cartesian(GCOORD_SPH); 
RR_X_90_CCW                = [ 1  0  0 ; ...
                               0  0 -1 ; ...
                               0  1  0 ];  % rotation matrix of 90° around X axis counterclockwise
RR_Z_180_CCW               = [-1  0  0 ; ...
                               0 -1  0 ; ...
                               0  0  1 ];   % rotation matrix of 180° around Z axis counterclockwise
r_cmb                      = MESH.r_cmb;  % CMB radius (r min)
r_surf                     = MESH.r_surf; % surface radius (r max)
r_tol                      = 1e-12 * (r_surf - r_cmb); % Tolerance to place points INSIDE the mesh (rather then exaclty ON the boundary)


%===========================================================================================================================================================================================================================================
% LOCATE ALL POINTS IN ELEMENTS INSIDE THE CONE, i.e., IN THE 90° X-AXIS ROTATED SPHERICAL FRAME
%===========================================================================================================================================================================================================================================
gX_PT                                    = spherical2cartesian(gTH_PT);
GCOORD_rot_X_90                          = RR_X_90_CCW * GCOORD;                 % rotate points of the mesh
GCOORD_SPH_rot_X_90                      = cartesian2spherical(GCOORD_rot_X_90); % transform to spherical coordinates
GCOORD_SPH_rot_X_90(3,MESH.PointID==301) = MESH.r_cmb;  % make sure that cmb nodes are exactly on the cmb (the radius of cmb nodes can shift slightly (~1e-8) from the actual cmb due to coordinate transformation)
GCOORD_SPH_rot_X_90(3,MESH.PointID==306) = MESH.r_surf; % make sure that surface nodes are exactly on the surface (the radius of surface nodes can shift slightly (~1e-8) from the actual surface due to coordinate transformation)
gX_PT_rot_X_90                           = RR_X_90_CCW * gX_PT;                  % rotate points to be located
gTH_PT_rot_X_90                          = cartesian2spherical(gX_PT_rot_X_90);  % transform to spherical coordinates
gTH_PT_rot_X_90(3,gTH_PT_rot_X_90(3,:) <= MESH.r_cmb + r_tol)  = MESH.r_cmb + r_tol;  % some cmb nodes may have been shifted outside the cmb by a tiny amount (~1e-12) due to the transformations (rotation + coordinate transformation) -> move them to the cmb 
gTH_PT_rot_X_90(3,gTH_PT_rot_X_90(3,:) >= MESH.r_surf - r_tol) = MESH.r_surf - r_tol; % some surface nodes may have been shifted outside the surface by a tiny amount (~1e-12) due to the transformations (rotation + coordinate transformation) -> move them to the surface 
[els_in_cone_rot_X_90,~,~] = ...
    tsearch2(GCOORD_SPH_rot_X_90(:,1:nnod),EL2NOD([1 2 4 3],:),gTH_PT_rot_X_90,WS_els_in_cone_rot_X_90,[]);
iloc                       = ismember(els_in_cone_rot_X_90,union(els_in_cone_iso,els_in_cone_no_iso));
els(iloc,1)                = els_in_cone_rot_X_90(iloc)';

% lightMag   = 0.8*[1 0.7 1]; % colour for the elements within the cone
% lightWhi   = [1 1 1];       % colour for elements outside the cone and not crossing phi = 2pi
% lightYell  = 0.8*[1 1 0];   % colour for elements outside the cone and crossing phi = 2pi
% lightCyan  = 0.8*[0 1 1];   % colour for the elements within the cone and isoparametric
% lightGreen = [0.5 0.8 0];   % colour for located nodes
% lightRed   = [0.8 0 0];     % colour for unlocated nodes
% % Points located
% plot_points_and_elements(GCOORD,GCOORD_SPH_rot_X_90,EL2NOD,gX_PT(:,iloc),gTH_PT_rot_X_90(:,iloc),els(els(:,1)~=0,1),45,lightMag,lightGreen,200)
% subplot(1,2,1)
% title('Inside-cone points located in elements inside the cone in Cartesian coordinates')
% subplot(1,2,2)
% title('Inside-cone points located in elements inside the cone in 90° rotated spherical coordinates')
% % Points not located
% plot_points_and_elements(GCOORD,GCOORD_SPH_rot_X_90,EL2NOD,gX_PT(:,~iloc),gTH_PT_rot_X_90(:,~iloc),els(els(:,1)~=0,1),45,lightMag,lightRed,300)
% subplot(1,2,1)
% title('Points still NOT located in Cartesian coordinates (MAGENTA ELS = elements inside the cone)')
% subplot(1,2,2)
% title('Points still NOT located in 90° rotated spherical coordinates (MAGENTA ELS = elements inside the cone)')

%===========================================================================================================================================================================================================================================
% LOCATE OUTSIDE-CONE POINTS IN ELEMENTS OUTSIDE THE CONE, i.e., IN THE ORIGINAL SPHERICAL FRAME
%===========================================================================================================================================================================================================================================
gTH_PT(3,gTH_PT(3,:) - MESH.r_cmb < 0)  = MESH.r_cmb;
gTH_PT(3,gTH_PT(3,:) - MESH.r_surf > 0) = MESH.r_surf;
gTH_PT_out_cone                         = gTH_PT(:,points_out_cone);  % select points to be located that are outside the cone
[els_out_cone,~,~]   = ...
    tsearch2(GCOORD_SPH(:,1:nnod),EL2NOD([1 2 4 3],:),gTH_PT_out_cone,WS_els_out_cone,[]);
iloc                 = ismember(els_out_cone,els_out_cone_no_cross_2pi);
iloc_out_cone        = points_out_cone(iloc);
els(iloc_out_cone,2) = els_out_cone(iloc)';
points_out_cone_left = points_out_cone(~iloc);

% % Points located
% gX_PT_out_cone = gX_PT(:,points_out_cone);
% plot_points_and_elements(GCOORD,GCOORD_SPH,EL2NOD,gX_PT_out_cone(:,iloc),gTH_PT_out_cone(:,iloc),els(els(:,2)~=0,2),45,lightWhi,lightGreen,201)
% subplot(1,2,1)
% title('Outside-cone points located in elements outside the cone and not crossing \phi = 2\pi in Cartesian coordinates')
% subplot(1,2,2)
% title('Outside-cone points located in elements outside the cone and not crossing \phi = 2\pi in spherical coordinates')
% % Points not located
% plot_points_and_elements(GCOORD,GCOORD_SPH,EL2NOD,gX_PT_out_cone(:,~iloc),gTH_PT_out_cone(:,~iloc),els(els(:,2)~=0,2),45,lightWhi,lightRed,301)
% subplot(1,2,1)
% title('Points still NOT located in Cartesian coordinates (WHITE ELS = elements outside the cone not crossing \phi = 2\pi)')
% subplot(1,2,2)
% title('Points still NOT located in spherical coordinates (WHITE ELS = elements outside the cone not crossing \phi = 2\pi)')

%===========================================================================================================================================================================================================================================
% LOCATE REMAINING OUTSIDE POINTS TO BE LOCATED IN ELEMENTS CROSSING PHI = 2PI, i.e., IN THE 180° Z-AXIS ROTATED SPHERICAL FRAME
%===========================================================================================================================================================================================================================================
GCOORD_rot_Z_180                          = RR_Z_180_CCW * GCOORD;                 % rotate points of the mesh
GCOORD_SPH_rot_Z_180                      = cartesian2spherical(GCOORD_rot_Z_180); % transform to spherical coordinates
GCOORD_SPH_rot_Z_180(3,MESH.PointID==301) = MESH.r_cmb;  % make sure that cmb nodes are exactly on the cmb (the radius of cmb nodes can shift slightly (~1e-8) from the actual cmb due to coordinate transformation)
GCOORD_SPH_rot_Z_180(3,MESH.PointID==306) = MESH.r_surf; % make sure that surface nodes are exactly on the surface (the radius of surface nodes can shift slightly (~1e-8) from the actual surface due to coordinate transformation)
gX_PT_rot_Z_180                           = RR_Z_180_CCW * gX_PT;                  % rotate points to be located
gTH_PT_rot_Z_180                          = cartesian2spherical(gX_PT_rot_Z_180);  % transform to spherical coordinates
gTH_PT_rot_Z_180(3,gTH_PT_rot_Z_180(3,:) <= MESH.r_cmb + r_tol)  = MESH.r_cmb + r_tol;  % some cmb nodes may have been shifted outside the cmb by a tiny amount (~1e-12) due to the transformations (rotation + coordinate transformation) -> move them to the cmb 
gTH_PT_rot_Z_180(3,gTH_PT_rot_Z_180(3,:) >= MESH.r_surf - r_tol) = MESH.r_surf - r_tol; % some surface nodes may have been shifted outside the surface by a tiny amount (~1e-12) due to the transformations (rotation + coordinate transformation) -> move them to the surface
[els_cross_2pi_rot_Z_180,~,~] = ...
    tsearch2(GCOORD_SPH_rot_Z_180(:,1:nnod),EL2NOD([1 2 4 3],:),gTH_PT_rot_Z_180(:,points_out_cone_left),WS_els_cross_2pi_rot_Z_180,[]);
iloc                  = ismember(els_cross_2pi_rot_Z_180,els_out_cone_cross_2pi);
iloc_out_cone         = points_out_cone_left(iloc);
els(iloc_out_cone,3)  = els_cross_2pi_rot_Z_180(iloc)';
points_out_cone_left2 = points_out_cone_left(~iloc);

% % Points located
% gX_PT_cross_2pi            = gX_PT(:,points_out_cone_left);
% gTH_PT_rot_Z_180_cross_2pi = gTH_PT_rot_Z_180(:,points_out_cone_left);
% plot_points_and_elements(GCOORD,GCOORD_SPH_rot_Z_180,EL2NOD,gX_PT_cross_2pi(:,iloc),gTH_PT_rot_Z_180_cross_2pi(:,iloc),els(els(:,3)~=0,3),45,lightYell,lightGreen,202)
% subplot(1,2,1)
% title('Outside-cone points located in elements outside the cone and crossing \phi = 2\pi in Cartesian coordinates')
% subplot(1,2,2)
% title('Outside-cone points located in elements outside the cone and crossing \phi = 2\pi in 180° rotated spherical coordinates')
% % Points not located
% plot_points_and_elements(GCOORD,GCOORD_SPH_rot_Z_180,EL2NOD,gX_PT_cross_2pi(:,~iloc),gTH_PT_rot_Z_180_cross_2pi(:,~iloc),els(els(:,3)~=0,3),45,lightYell,lightRed,302)
% subplot(1,2,1)
% title('Points still NOT located in Cartesian coordinates (YELLOW ELS = elements outside the cone crossing \phi = 2\pi)')
% subplot(1,2,2)
% title('Points still NOT located in 180° rotated spherical coordinates (YELLOW ELS = elements outside the cone crossing \phi = 2\pi)')

%===========================================================================================================================================================================================================================================
% LOCATE REMAINING OUTSIDE POINTS IN ELEMENTS INSIDE THE CONE, i.e., IN THE 90° X-AXIS ROTATED SPHERICAL FRAME
%===========================================================================================================================================================================================================================================
gTH_PT_out_cone_rot_X_90   = gTH_PT_rot_X_90(:,points_out_cone_left2); % select points to be located that are outside the cone in the rotated sph frame
[els_in_cone_rot_X_90,~,~] = ...
    tsearch2(GCOORD_SPH_rot_X_90(:,1:nnod),EL2NOD([1 2 4 3],:),gTH_PT_out_cone_rot_X_90,WS_els_in_cone_rot_X_90,[]);
iloc                       = els_in_cone_rot_X_90 > 0;
iloc2                      = points_out_cone_left2(iloc);
iloc3                      = els(iloc2,1)==0;
els(iloc2(iloc3),4)        = els_in_cone_rot_X_90(iloc3)';
points_out_cone_left3      = points_out_cone_left2(~iloc);

% % Points located
% gX_PT_out_cone_rot_X_90 = gX_PT(:,points_out_cone_left2);
% plot_points_and_elements(GCOORD,GCOORD_SPH_rot_X_90,EL2NOD,gX_PT_out_cone_rot_X_90(:,iloc),gTH_PT_out_cone_rot_X_90(:,iloc),els(els(:,4)~=0,4),45,lightCyan,lightGreen,203)
% subplot(1,2,1)
% title('Outside-cone points located in isoparametric elements in Cartesian coordinates')
% subplot(1,2,2)
% title('Outside-cone points located in isoparametric elements in 90° rotated spherical coordinates')
% % Points not located
% plot_points_and_elements(GCOORD,GCOORD_SPH_rot_X_90,EL2NOD,gX_PT_out_cone_rot_X_90(:,~iloc),gTH_PT_out_cone_rot_X_90(:,~iloc),els(els(:,4)~=0,4),45,lightCyan,lightRed,303)
% subplot(1,2,1)
% title('Points still NOT located in Cartesian coordinates (CYAN ELS = isoparametric elements inside the cone)')
% subplot(1,2,2)
% title('Points still NOT located in 90° rotated spherical coordinates (CYAN ELS = isoparametric elements inside the cone)')

if ~isempty(points_out_cone_left3)
    %=======================================================================================================================================================================================================================================
    % LOCATE REMAINING OUTSIDE POINTS IN ELEMENTS INSIDE THE CONE IN THE ORIGINAL SPHERICAL FRAME (special case in which the element is inside the cone with just one vertex and it is not isoparametric) 
    %=======================================================================================================================================================================================================================================
    % The vertical cross-section would be:
    %
    %   \                                    /
    %    \                                  /
    %     \                                /
    %      \         inside  cone         /      outside  cone
    %       \                            /      
    %        \                          /
    %         \                        /
    %          \                 - - -/- - - - - 
    %           \                \   /         /   <-- This tetrahedron is considered inside the cone since it has one of the vertices inside the cone.
    %            \                \ / x       /        However, its edges are straight, since the remaining vertices are outside the cone, so it is not an isoparametric element 
    %             \                /         /         (remember that isoparametric elements have at least one bar inside the cone, i.e., 2 vertices inside the cone).  
    %              \              / \       /          So, this tetrahedron is 'inside' the cone but its edges are straight in the original spherical frame 
    %               \            /   \     /
    %                \          /     \   /
    %                 \        /       \ /
    %                  \      /
    
    gTH_PT_out_cone2      = gTH_PT(:,points_out_cone_left3); % select points to be located that are outside the cone in the original sph frame
    [els_in_cone,~,~]     = tsearch2(GCOORD_SPH(:,1:nnod),EL2NOD([1 2 4 3],:),gTH_PT_out_cone2,[],[]);
    iloc                  = els_in_cone > 0;
    iloc2                 = points_out_cone_left3(iloc);
    iloc3                 = els(iloc2,1)==0;
    els                   = [els zeros(size(gTH_PT,2),1)];
    els(iloc2(iloc3),5)   = els_in_cone(iloc3)';
    points_out_cone_left4 = points_out_cone_left3(~iloc);
    
%     % Points located
%     gX_PT_out_cone2 = gX_PT(:,points_out_cone_left3);
%     plot_points_and_elements(GCOORD,GCOORD_SPH,EL2NOD,gX_PT_out_cone2(:,iloc),gTH_PT_out_cone2(:,iloc),els(els(:,5)~=0,5),45,lightCyan,lightGreen,204)
%     subplot(1,2,1)
%     title('Outside-cone points located in isoparametric elements in Cartesian coordinates')
%     subplot(1,2,2)
%     title('Outside-cone points located in isoparametric elements in the original spherical coordinates')
%     % Points not located
%     plot_points_and_elements(GCOORD,GCOORD_SPH,EL2NOD,gX_PT_out_cone2(:,~iloc),gTH_PT_out_cone2(:,~iloc),els(els(:,5)~=0,5),45,lightCyan,lightRed,304)
%     subplot(1,2,1)
%     title('Points still NOT located in Cartesian coordinates (CYAN ELS = isoparametric elements inside the cone)')
%     subplot(1,2,2)
%     title('Points still NOT located in the original spherical coordinates (CYAN ELS = isoparametric elements inside the cone)')
    
    if ~isempty(points_out_cone_left4)
        
        %=======================================================================================================================================================================================================================================
        % LOCATE REMAINING OUTSIDE POINTS IN ELEMENTS INSIDE THE CONE IN THE 180° Z-AXIS ROTATED SPHERICAL FRAME (special case in which the element is inside the cone with just one vertex and it is not isoparametric) 
        %=======================================================================================================================================================================================================================================
        % The vertical cross-section would be:
        %
        %   \                                    /
        %    \                                  /
        %     \                                /
        %      \         inside  cone         /      outside  cone
        %       \                            /      
        %        \                          /
        %         \                        /
        %          \                 - - -/- - - - - 
        %           \                \   /         /   <-- This tetrahedron is considered inside the cone since it has one of the vertices inside the cone.
        %            \                \ / x       /        However, its edges are straight, since the remaining vertices are outside the cone, so it is not an isoparametric element 
        %             \                /         /         (remember that isoparametric elements have at least one bar inside the cone, i.e., 2 vertices inside the cone).  
        %              \              / \       /          So, this tetrahedron is 'inside' the cone but its edges are straight. Furthermore, since the tetraheron is crossing 2pi 
        %               \            /   \     /           we need to search in the 180 Z-axis rotated spherical frame 
        %                \          /     \   /
        %                 \        /       \ /
        %                  \      /
        
        gTH_PT_out_cone3      = gTH_PT_rot_Z_180(:,points_out_cone_left4);
        [els_in_cone2,~,~]    = tsearch2(GCOORD_SPH_rot_Z_180(:,1:nnod),EL2NOD([1 2 4 3],:),gTH_PT_out_cone3,[],[]);
        iloc                  = els_in_cone2 > 0;
        iloc2                 = points_out_cone_left4(iloc);
        iloc3                 = els(iloc2,1)==0;
        els                   = [els zeros(size(gTH_PT,2),1)];
        els(iloc2(iloc3),5)   = els_in_cone2(iloc3)';
        points_out_cone_left5 = points_out_cone_left4(~iloc);
        
%         % Points located
%         gX_PT_out_cone3 = gX_PT(:,points_out_cone_left4);
%         plot_points_and_elements(GCOORD,GCOORD_SPH,EL2NOD,gX_PT_out_cone3(:,iloc),gTH_PT_out_cone3(:,iloc),els(els(:,5)~=0,5),45,lightCyan,lightGreen,205)
%         subplot(1,2,1)
%         title('Outside-cone points located in isoparametric elements in Cartesian coordinates')
%         subplot(1,2,2)
%         title('Outside-cone points located in isoparametric elements in the 180 Z-axis rotated spherical coordinates')
%         % Points not located
%         plot_points_and_elements(GCOORD,GCOORD_SPH,EL2NOD,gX_PT_out_cone3(:,~iloc),gTH_PT_out_cone3(:,~iloc),els(els(:,5)~=0,5),45,lightCyan,lightRed,305)
%         subplot(1,2,1)
%         title('Points still NOT located in Cartesian coordinates (CYAN ELS = isoparametric elements inside the cone)')
%         subplot(1,2,2)
%         title('Points still NOT located in the 180 Z-axis rotated spherical coordinates (CYAN ELS = isoparametric elements inside the cone)')
        
        if ~isempty(points_out_cone_left5)
            error('no good');
        end
    end
end

if min(sum(els,2))==0
    %=======================================================================================================================================================================================================================================
    % LOCATE REMAINING INSIDE POINTS IN ELEMENTS INSIDE THE CONE IN THE ORIGINAL SPHERICAL FRAME (special case in which the element is inside the cone with just one vertex and it is not isoparametric) 
    %=======================================================================================================================================================================================================================================
    % The vertical cross-section would be:
    %
    %   \                                    /
    %    \                                  /
    %     \                                /
    %      \         inside  cone         /      outside  cone
    %       \                            /      
    %        \                          /
    %         \                        /
    %          \                 - - -/- - - - - 
    %           \                \ x /         /   <-- This tetrahedron is considered inside the cone since it has one of the vertices inside the cone.
    %            \                \ /         /        However, its edges are straight, since the remaining vertices are outside the cone, so it is not an isoparametric element 
    %             \                /         /         (remember that isoparametric elements have at least one bar inside the cone, i.e., 2 vertices inside the cone).  
    %              \              / \       /          So, this tetrahedron is 'inside' the cone but its edges are straight in the original spherical frame 
    %               \            /   \     /
    %                \          /     \   /
    %                 \        /       \ /
    %                  \      /
    
    points_in_cone_left  = find(sum(els,2)==0);
    gTH_PT_in_cone       = gTH_PT(:,points_in_cone_left); % select points to be located that are inside the cone in the original sph frame
    [els_in_cone2,~,~]   = tsearch2(GCOORD_SPH(:,1:nnod),EL2NOD([1 2 4 3],:),gTH_PT_in_cone,[],[]);
    iloc4                = ismember(els_in_cone2,union(els_in_cone_iso,els_in_cone_no_iso));
    iloc5                = points_in_cone_left(iloc4);
    iloc6                = els(iloc5,1)==0;
    els(iloc5(iloc6),2)  = els_in_cone2(iloc4)';
    points_in_cone_left2 = points_in_cone_left(~iloc4);
    
    if ~isempty(points_in_cone_left2)
        %=======================================================================================================================================================================================================================================
        % LOCATE REMAINING INSIDE POINTS IN ELEMENTS INSIDE THE CONE IN THE 180� Z-AXIS ROTATED SPHERICAL FRAME (special case in which the element is inside the cone with just one vertex and it is not isoparametric)
        %=======================================================================================================================================================================================================================================
        % The vertical cross-section would be:
        %
        %   \                                    /
        %    \                                  /
        %     \                                /
        %      \         inside  cone         /      outside  cone
        %       \                            /
        %        \                          /
        %         \                        /
        %          \                 - - -/- - - - -
        %           \                \ x /         /   <-- This tetrahedron is considered inside the cone since it has one of the vertices inside the cone.
        %            \                \ /         /        However, its edges are straight, since the remaining vertices are outside the cone, so it is not an isoparametric element
        %             \                /         /         (remember that isoparametric elements have at least one bar inside the cone, i.e., 2 vertices inside the cone).
        %              \              / \       /          So, this tetrahedron is 'inside' the cone but its edges are straight. Furthermore, since the tetraheron is crossing 2pi 
        %               \            /   \     /           we need to search in the 180 Z-axis rotated spherical frame 
        %                \          /     \   /
        %                 \        /       \ /
        %                  \      /
        
        gTH_PT_in_cone2      = gTH_PT_rot_Z_180(:,points_in_cone_left2);
        [els_in_cone3,~,~]   = tsearch2(GCOORD_SPH_rot_Z_180(:,1:nnod),EL2NOD([1 2 4 3],:),gTH_PT_in_cone2,[],[]);
        iloc7                = ismember(els_in_cone3,union(els_out_cone_cross_2pi,els_out_cone_no_cross_2pi));
        iloc8                = points_in_cone_left2(iloc7);
        iloc9                = els(iloc8,1)==0;
        els(iloc8(iloc9),2)  = els_in_cone3(iloc7)';
        points_in_cone_left3 = points_in_cone_left2(~iloc7);
        
        if ~isempty(points_in_cone_left3)
            % It might happen that a point which is inside the cone was not well located in the 90 X-axis rotated spherical frame.
            % This is because this point, in the 90 X-axis rotated spherical frame, might be found (using tsearch2) in two different elements: 
            % - one of them is the right element
            % - the other one is an element that is crossing 2pi in this  90 X-axis rotated spherical frame (some vertices of the element 
            %   are close to 0 and some other are close to 2pi making (when plotting in theta, phi, r) a stretched element which is not the good one
            % If this situation happens, tsearch2 will assign the number of the first element where the point is located, even if it is not the right element 
            % Then, when using iloc = ismember(els_in_cone_rot_X_90,union(els_in_cone_iso,els_in_cone_no_iso)) the point is not located inside els_in_cone_iso
            % or els_in_cone_no_iso because tsearch assigned previously an element that was crossing 2pi
            % To solve this issue, use tsarch2 to find the element (even if is not the right one -> it will corrected using the routine locate_points_in_neighboring_element 
            [els_in_cone_rot_X_90_v2,~,~] = ...
                tsearch2(GCOORD_SPH_rot_X_90(:,1:nnod),EL2NOD([1 2 4 3],:),gTH_PT_rot_X_90(:,points_in_cone_left3),WS_els_in_cone_rot_X_90,[]);
            iloc10                = els_in_cone_rot_X_90_v2 > 0;
            iloc11                = points_in_cone_left3(iloc10);
            iloc12                = els(iloc11,1)==0;
            els(iloc11(iloc12),1) = els_in_cone_rot_X_90_v2(iloc10)';
            points_in_cone_left4  = points_in_cone_left3(~iloc10);
            if ~isempty(points_in_cone_left4)
                error('no good');
            end
        end
    end
end

%===========================================================================================================================================================================================================================================
% OUTPUT
%===========================================================================================================================================================================================================================================
ipt_problem = find(sum(els>0,2)>1);
if ~isempty(ipt_problem)
    for j=1:length(ipt_problem)
        ind          = find(els(ipt_problem(j),:)>0);
        els_possible = els(ipt_problem(j),ind);
        lc_possible  = local_coords_3d_sph(GCOORD_SPH,EL2NOD,els_possible,...
            repmat(gTH_PT(:,ipt_problem(j)),1,length(els_possible)),MESH);
        lc_possible(4,:) = 1 - sum(lc_possible,1);
        ilc_wrng         = min(lc_possible,[],1)<0 | max(lc_possible,[],1)>1;
        if all(~ilc_wrng) && sum(ismember(els_possible,MESH.els_in_cone_iso)) == 1
            ilc_wrng = ismember(els_possible,MESH.els_in_cone_iso);
        end
        if all(~ilc_wrng) && all(min(lc_possible) == zeros(1,size(ind,2))) && all(max(lc_possible) == ones(1,size(ind,2)))
            % it means the the point is in a shared vertex of two (or three) elements, pick just one 
            i_w      = [0 1];
            ilc_wrng = i_w == 1;
        end
        if all(~ilc_wrng) && sum(ismember(els_possible,MESH.els_in_cone_iso)) == 0
            error('no good');
        end
        if all(ilc_wrng)
            error('no good');
        end
        els(ipt_problem(j),ind(ilc_wrng)) = 0;
    end
%    
end
if min(sum(els,2))==0
    error('no good either.');
end
els = uint32(max(els,[],2));

%===========================================================================================================================================================================================================================================
% CHECK IF THERE IS STILL ANY POINT NOT LOCATED
%===========================================================================================================================================================================================================================================
iloc  = ~isnan(els) & els>0;
ilost = find(~iloc);
if any(els(ilost)==0)
    figure(777)
    clf
    scatter3(gTH_PT(1,ilost),gTH_PT(2,ilost),gTH_PT(3,ilost)/1e3,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
    axis equal
    axis([0 pi 0 2*pi 3 6.5])
    view(142.5,30)
    xlabel('\theta (rad)')
    ylabel('\phi (rad)')
    zlabel('r (10^3 km)')
    error('some nodes have not been located');
end

end % END OF FUNCTION tsearch2_sph

% ##########################################################################################################################################################################################################
%                                           SUB-FUNCTIONS
% ##########################################################################################################################################################################################################

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