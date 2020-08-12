function [GCOORD_cubic,EL2NOD_cubic] = cubic_els_in_cone_iso(GCOORD4,GCOORD,GCOORD4_SPH,EL2NOD_c,els_in_cone_iso,SETTINGS)
% Usage: [GCOORD_cubic,EL2NOD_cubic] = cubic_els_in_cone_iso(GCOORD4,GCOORD,GCOORD4_SPH,EL2NOD_c,els_in_cone_iso,SETTINGS)
%
% Purpose: 
%   Compute 20 nodes and connectivity for isoparametric elements. The
%   mid-face points are computed using the information of the 9 nodes (3
%   vertices + 6 edge points) that define each face.
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
%   GCOORD4         : [matrix]    : Cartesian coordinates of all nodes 
%                                   (4 nodel) (3 x nnod4)
%   GCOORD          : [matrix]    : Cartesian coordinates of all nodes 
%                                   (10 nodel) (3 x nnod10)
%   GCOORD4_SPH     : [matrix]    : spherical coordinates of all nodes 
%                                   (4 nodel) (3 x nnod4)
%   EL2NOD_c        : [matrix]    : connectivity matrix (4 x nel)
%   els_in_cone_iso : [vector]    : isoparametric elements
%   SETTINGS        : [structure] : model parameters
%
% Output:
%   GCOORD_cubic    : [matrix]    : Cartesian coordinates of all edges 
%                                   (1/3 and 2/3) and face nodes  
%   EL2NOD_cubic    : [matrix]    : connectivity matrix for isoparametric
%                                   elements (20 x nels_in_cone_iso)
%
% JMT Aug 2017
%

%=========================================================================================================================
% COMPUTE NODES AT 1/3 and 2/3 IN THE EDGES IN THE ORIGINAL FRAME
%=========================================================================================================================
% Size of coarse mesh
EL2NOD4         = EL2NOD_c(1:4,els_in_cone_iso); % EL2NOD4 for isoparametric elements
nel             = size(EL2NOD4,2);             % number of elements
nnod4           = size(GCOORD4_SPH,2);         % number of nodes
% Create a pointer edges defined by their end-nodes (e.g. bar 1 has end-nodes [1 2], bar 2 has [2 3],...)
edges           = [ 1  2 ... % 1st edge of parent
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
thBarEnds       = reshape(GCOORD4_SPH(1,EDGE2NOD'),2,nedge);
phBarEnds       = reshape(GCOORD4_SPH(2,EDGE2NOD'),2,nedge);
rBarEnds        = reshape(GCOORD4_SPH(3,EDGE2NOD'),2,nedge);

% Check if phi angles (for bar nodes) are separated an angular distance bigger than pi 
phBarEnds       = check_ang_dist_phi(GCOORD4_SPH,EDGE2NOD,phBarEnds);

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

% storage for cubic order connectivity matrix
EL2NOD20         = zeros(16,nel,'uint32'); % must be 16, face and bubble nodes added later
EL2NOD20(1: 4,:) = EL2NOD4;
EL2NOD20(5:16,:) = nnod4 + reshape([2*EL2EDGE(:)'-1; 2*EL2EDGE(:)'],12,[]);

%======================================================================================================================
% Now correct local node numbering for edge nodes (ouside the cone crossing 2pi and outside the cone not crossing 2pi) 
%======================================================================================================================
% Select those edges outside the cone boundary (both edge ends are outside of the cone boundary)
TH2EDGES               = reshape(GCOORD4_SPH(1,EDGE2NOD),[],2);
theta_bars_out_cone    = ...
    sum(TH2EDGES > SETTINGS.theta_cone*pi/180 & TH2EDGES < (180-SETTINGS.theta_cone)*pi/180,2); % vector for edges having both ends outside the cone boundary
                                         % 0 --> both ends bar are inside the cone boundary
                                         % 1 --> one end bar is inside the cone boundary and the other one is outside the cone boundary
                                         % 2 --> both ends bar are outside the cone boundary
EDGE2NOD_out_cone                     = EDGE2NOD(theta_bars_out_cone == 2,:);
ang_dist                              = abs(GCOORD4_SPH(2,EDGE2NOD_out_cone(:,1))-GCOORD4_SPH(2,EDGE2NOD_out_cone(:,2)));
ind                                   = ang_dist > pi; % indices for those bars having the phi angles (of the ends) separated more than pi
EDGE2NOD_out_cone_cross_2pi           = EDGE2NOD_out_cone(ind,:); % edges ouside the cone boundary and crossing 2pi 
EDGE2NOD_out_cone_no_cross_2pi        = EDGE2NOD_out_cone;
EDGE2NOD_out_cone_no_cross_2pi(ind,:) = []; % edges ouside the cone boundary and not crossing 2pi

% Coordinates for vertices (4) and nodes on the edges (12) in original frame
nnod16                              = nnod4 + 2*nedge;
GCOORD16_SPH_temp                   = zeros(3,nnod16);
GCOORD16_SPH_temp(:,1:nnod4)        = GCOORD4_SPH;
GCOORD16_SPH_temp(1,nnod4+1:nnod16) = thBarThirds;
GCOORD16_SPH_temp(2,nnod4+1:nnod16) = phBarThirds;
GCOORD16_SPH_temp(3,nnod4+1:nnod16) = rBarThirds;
GCOORD16_temp                       = spherical2cartesian(GCOORD16_SPH_temp); % convert to Cartesian coordinates

% Coordinates for vertices (4) and nodes on the edges (12) in the 180° rotated around Z axis frame
RR_Z_180_CCW      = [-1  0  0 ; 0 -1  0 ; 0  0  1 ];    % 180° counterclockwise around the Z axis rotation matrix
GCOORD_rot180     = RR_Z_180_CCW * GCOORD16_temp;       % rotate coordinates
GCOORD_SPH_rot180 = cartesian2spherical(GCOORD_rot180); % convert to spherical coordinates

% local coordinates of all edge nodes (nodes 5:16)
nvertx = 4;
r      = [2/3 1/3  0   0   0   0  1/3 2/3 2/3 1/3  0   0 ];
s      = [1/3 2/3 2/3 1/3  0   0   0   0   0   0  1/3 2/3];
t      = [ 0   0  1/3 2/3 2/3 1/3  0   0  1/3 2/3  0   0 ];
N      = sf_dsf_tet([r;s;t],nvertx,'matrix');

% Compute (theta,phi,r) given by the mesh coordinates and (theta,phi,r) given by the shape functions in original frame 
thn_vertx         = reshape(GCOORD16_SPH_temp(1,EL2NOD20(1:nvertx,:)),nvertx,[]);
thn_check         = N' * thn_vertx;
thn_mesh          = reshape(GCOORD16_SPH_temp(1,EL2NOD20(5:16,:)),12,[]);
phn_vertx         = reshape(GCOORD16_SPH_temp(2,EL2NOD20(1:nvertx,:)),nvertx,[]);
phn_check         = N' * phn_vertx;
phn_mesh          = reshape(GCOORD16_SPH_temp(2,EL2NOD20(5:16,:)),12,[]);
rn_vertx          = reshape(GCOORD16_SPH_temp(3,EL2NOD20(1:nvertx,:)),nvertx,[]);
rn_check          = N' * rn_vertx;
rn_mesh           = reshape(GCOORD16_SPH_temp(3,EL2NOD20(5:16,:)),12,[]);

% Compute (theta,phi,r) given by the mesh coordinates and (theta,phi,r) given by the shape functions in the 180° rotated around Z axis frame 
thn_vertx_rot180  = reshape(GCOORD_SPH_rot180(1,EL2NOD20(1:nvertx,:)),nvertx,[]);
thn_check_rot180  = N' * thn_vertx_rot180;
thn_mesh_rot180   = reshape(GCOORD_SPH_rot180(1,EL2NOD20(5:16,:)),12,[]);
phn_vertx_rot180  = reshape(GCOORD_SPH_rot180(2,EL2NOD20(1:nvertx,:)),nvertx,[]);
phn_check_rot180  = N' * phn_vertx_rot180;
phn_mesh_rot180   = reshape(GCOORD_SPH_rot180(2,EL2NOD20(5:16,:)),12,[]);
rn_vertx_rot180   = reshape(GCOORD_SPH_rot180(3,EL2NOD20(1:nvertx,:)),nvertx,[]);
rn_check_rot180   = N' * rn_vertx_rot180;
rn_mesh_rot180    = reshape(GCOORD_SPH_rot180(3,EL2NOD20(5:16,:)),12,[]);

% Compare coordinates for each element
for i=1:nel
    this_el                      = EL2NOD4(:,i);
    edges_this_el                = reshape( this_el(edges',:),2,[] )'; % edges of this element
    i_out_cone_no_cross_2pi      = ismember(edges_this_el,EDGE2NOD_out_cone_no_cross_2pi,'rows');
    i_out_cone_no_cross_2pi_flip = ismember(fliplr(edges_this_el),EDGE2NOD_out_cone_no_cross_2pi,'rows');
    i_out_cone_no_cross_2pi      = (i_out_cone_no_cross_2pi + i_out_cone_no_cross_2pi_flip) > 0;
    i_out_cone_cross_2pi         = ismember(edges_this_el,EDGE2NOD_out_cone_cross_2pi,'rows');
    i_out_cone_cross_2pi_flip    = ismember(fliplr(edges_this_el),EDGE2NOD_out_cone_cross_2pi,'rows');
    i_out_cone_cross_2pi         = (i_out_cone_cross_2pi + i_out_cone_cross_2pi_flip) > 0;
    
%     if i == 115 % for debugging
%         
%         sph_vertx  = [thn_vertx(:,i) phn_vertx(:,i) rn_vertx(:,i)]';
%         cart_vertx = spherical2cartesian(sph_vertx);
%         sph_check  = [thn_check(:,i) phn_check(:,i) rn_check(:,i)]';
%         cart_check = [cart_vertx spherical2cartesian(sph_check)];
%         sph_mesh   = [thn_mesh(:,i) phn_mesh(:,i) rn_mesh(:,i)]';
%         cart_mesh  = [cart_vertx spherical2cartesian(sph_mesh)];
%         
%         figure(1); clf
%         lightMag    = 0.80*[1 0.7 1]; % colour for the cone boundary
%         lightGrey   = 0.90*[1 1 1];   % colour for the shell boundaries
%         lightWhi    = [1 1 1];       % colour for elements
%         th          = SETTINGS.theta_cone*pi/180;
%         h           = 6371;
%         r           = h*tan(th);
%         [R,A]       = meshgrid(linspace(0,r,2),linspace(0,2*pi,21));
%         X           = R .* cos(A);
%         Y           = R .* sin(A);
%         Z           = R/tan(th);
%         surface(X,Y,Z,'FaceColor','none','EdgeColor',lightMag)
%         hold on
%         surface(X,Y,-Z,'FaceColor','none','EdgeColor',lightMag)
%         % Plot the shell
%         [x_sph,y_sph,z_sph] = sphere(20);
%         x_sph               = x_sph*3471;
%         y_sph               = y_sph*3471;
%         z_sph               = z_sph*3471;
%         axis equal
%         surface(x_sph,y_sph,z_sph,'FaceColor','none','EdgeColor',lightGrey)
%         [x_sph,y_sph,z_sph] = sphere(30);
%         x_sph               = x_sph*6371;
%         y_sph               = y_sph*6371;
%         z_sph               = z_sph*6371;
%         surface(x_sph,y_sph,z_sph,'FaceColor','none','EdgeColor',lightGrey)
%         axis([-6371 6371 -6371 6371 -6371 6371])
%         view(142.5,30)
%         grid on
%         hold on
%         tetramesh(EL2NOD4(:,115)',GCOORD','FaceColor',lightWhi,'FaceAlpha',0.3)
%         color_local_vertex_1 = [1 0 0]; % red
%         color_local_vertex_2 = [1 1 0]; % yellow
%         color_local_vertex_3 = [0 1 0]; % green
%         color_local_vertex_4 = [0 0 1]; % blue
%         scatter3(cart_check(1,1), ...
%                  cart_check(2,1), ...
%                  cart_check(3,1),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_1)
%         scatter3(cart_check(1,2), ...
%                  cart_check(2,2), ...
%                  cart_check(3,2),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_2)
%         scatter3(cart_check(1,3), ...
%                  cart_check(2,3), ...
%                  cart_check(3,3),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_3)
%         scatter3(cart_check(1,4), ...
%                  cart_check(2,4), ...
%                  cart_check(3,4),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_4)
%         
%         scatter3(cart_check(1,5), ...
%                  cart_check(2,5), ...
%                  cart_check(3,5),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_1)
%         scatter3(cart_check(1,6), ...
%                  cart_check(2,6), ...
%                  cart_check(3,6),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_2)
%         scatter3(cart_check(1,7), ...
%                  cart_check(2,7), ...
%                  cart_check(3,7),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_2)
%         scatter3(cart_check(1,8), ...
%                  cart_check(2,8), ...
%                  cart_check(3,8),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_3)
%         scatter3(cart_check(1,9), ...
%                  cart_check(2,9), ...
%                  cart_check(3,9),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_3)
%         scatter3(cart_check(1,10), ...
%                  cart_check(2,10), ...
%                  cart_check(3,10),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_4)
%         scatter3(cart_check(1,11), ...
%                  cart_check(2,11), ...
%                  cart_check(3,11),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_4)
%         scatter3(cart_check(1,12), ...
%                  cart_check(2,12), ...
%                  cart_check(3,12),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_1)
%         scatter3(cart_check(1,13), ...
%                  cart_check(2,13), ...
%                  cart_check(3,13),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_1)
%         scatter3(cart_check(1,14), ...
%                  cart_check(2,14), ...
%                  cart_check(3,14),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_3)
%         scatter3(cart_check(1,15), ...
%                  cart_check(2,15), ...
%                  cart_check(3,15),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_4)
%         scatter3(cart_check(1,16), ...
%                  cart_check(2,16), ...
%                  cart_check(3,16),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_2)
%         title('element check (Shape Functions)')
%         
%         
%         figure(2); clf
%         surface(X,Y,Z,'FaceColor','none','EdgeColor',lightMag)
%         hold on
%         surface(X,Y,-Z,'FaceColor','none','EdgeColor',lightMag)
%         % Plot the shell
%         [x_sph,y_sph,z_sph] = sphere(20);
%         x_sph               = x_sph*3471;
%         y_sph               = y_sph*3471;
%         z_sph               = z_sph*3471;
%         axis equal
%         surface(x_sph,y_sph,z_sph,'FaceColor','none','EdgeColor',lightGrey)
%         [x_sph,y_sph,z_sph] = sphere(30);
%         x_sph               = x_sph*6371;
%         y_sph               = y_sph*6371;
%         z_sph               = z_sph*6371;
%         surface(x_sph,y_sph,z_sph,'FaceColor','none','EdgeColor',lightGrey)
%         axis([-6371 6371 -6371 6371 -6371 6371])
%         view(142.5,30)
%         grid on
%         hold on
%         tetramesh(EL2NOD4(:,115)',GCOORD','FaceColor',lightWhi,'FaceAlpha',0.3)
%         scatter3(cart_mesh(1,1), ...
%                  cart_mesh(2,1), ...
%                  cart_mesh(3,1),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_1)
%         scatter3(cart_mesh(1,2), ...
%                  cart_mesh(2,2), ...
%                  cart_mesh(3,2),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_2)
%         scatter3(cart_mesh(1,3), ...
%                  cart_mesh(2,3), ...
%                  cart_mesh(3,3),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_3)
%         scatter3(cart_mesh(1,4), ...
%                  cart_mesh(2,4), ...
%                  cart_mesh(3,4),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_4)
%         
%         scatter3(cart_mesh(1,5), ...
%                  cart_mesh(2,5), ...
%                  cart_mesh(3,5),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_1)
%         scatter3(cart_mesh(1,6), ...
%                  cart_mesh(2,6), ...
%                  cart_mesh(3,6),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_2)
%         scatter3(cart_mesh(1,7), ...
%                  cart_mesh(2,7), ...
%                  cart_mesh(3,7),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_2)
%         scatter3(cart_mesh(1,8), ...
%                  cart_mesh(2,8), ...
%                  cart_mesh(3,8),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_3)
%         scatter3(cart_mesh(1,9), ...
%                  cart_mesh(2,9), ...
%                  cart_mesh(3,9),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_3)
%         scatter3(cart_mesh(1,10), ...
%                  cart_mesh(2,10), ...
%                  cart_mesh(3,10),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_4)
%         scatter3(cart_mesh(1,11), ...
%                  cart_mesh(2,11), ...
%                  cart_mesh(3,11),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_4)
%         scatter3(cart_mesh(1,12), ...
%                  cart_mesh(2,12), ...
%                  cart_mesh(3,12),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_1)
%         scatter3(cart_mesh(1,13), ...
%                  cart_mesh(2,13), ...
%                  cart_mesh(3,13),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_1)
%         scatter3(cart_mesh(1,14), ...
%                  cart_mesh(2,14), ...
%                  cart_mesh(3,14),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_3)
%         scatter3(cart_mesh(1,15), ...
%                  cart_mesh(2,15), ...
%                  cart_mesh(3,15),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_4)
%         scatter3(cart_mesh(1,16), ...
%                  cart_mesh(2,16), ...
%                  cart_mesh(3,16),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_2)
%         title('element mesh)')
%     
%     end
    
    if any(i_out_cone_no_cross_2pi)
        % Check local node numbering in edges ouside the cone and not crossing 2pi
        edges_out_cone_no_cross_2pi       = find(i_out_cone_no_cross_2pi > 0);
        nodes_edges_out_cone_no_cross_2pi = [2*edges_out_cone_no_cross_2pi-1 2*edges_out_cone_no_cross_2pi]'; 
        nodes_edges_out_cone_no_cross_2pi = nodes_edges_out_cone_no_cross_2pi(:); % it needs to be a colum vector
        
        diffth        = abs(thn_check(nodes_edges_out_cone_no_cross_2pi,i) - thn_mesh(nodes_edges_out_cone_no_cross_2pi,i)) > 1e-12;
        diffph        = abs(phn_check(nodes_edges_out_cone_no_cross_2pi,i) - phn_mesh(nodes_edges_out_cone_no_cross_2pi,i)) > 1e-12;
        diffr         = abs(rn_check(nodes_edges_out_cone_no_cross_2pi,i)  - rn_mesh(nodes_edges_out_cone_no_cross_2pi,i))  > 1e-08;
        
        % Correct local node numbering in edges ouside the cone and not crossing 2pi
        for inod = 1:2:size(nodes_edges_out_cone_no_cross_2pi,1)
            if (diffth(inod) || diffph(inod) || diffr(inod))
                EL2NOD20([4+nodes_edges_out_cone_no_cross_2pi(inod) 4+nodes_edges_out_cone_no_cross_2pi(inod+1)],i) = ...
                    EL2NOD20([4+nodes_edges_out_cone_no_cross_2pi(inod+1) 4+nodes_edges_out_cone_no_cross_2pi(inod)],i);
            end
        end
    end
    
%     if i == 115 % for debugging
%         thn_mesh_corrected  = reshape(GCOORD16_SPH_temp(1,EL2NOD20(5:16,:)),12,[]);
%         phn_mesh_corrected  = reshape(GCOORD16_SPH_temp(2,EL2NOD20(5:16,:)),12,[]);
%         rn_mesh_corrected   = reshape(GCOORD16_SPH_temp(3,EL2NOD20(5:16,:)),12,[]);
%         sph_mesh_corrected  = [thn_mesh_corrected(:,i) phn_mesh_corrected(:,i) rn_mesh_corrected(:,i)]';
%         cart_mesh_corrected = [cart_vertx spherical2cartesian(sph_mesh_corrected)];
%         
%         figure(3); clf
%         surface(X,Y,Z,'FaceColor','none','EdgeColor',lightMag)
%         hold on
%         surface(X,Y,-Z,'FaceColor','none','EdgeColor',lightMag)
%         % Plot the shell
%         [x_sph,y_sph,z_sph] = sphere(20);
%         x_sph               = x_sph*3471;
%         y_sph               = y_sph*3471;
%         z_sph               = z_sph*3471;
%         axis equal
%         surface(x_sph,y_sph,z_sph,'FaceColor','none','EdgeColor',lightGrey)
%         [x_sph,y_sph,z_sph] = sphere(30);
%         x_sph               = x_sph*6371;
%         y_sph               = y_sph*6371;
%         z_sph               = z_sph*6371;
%         surface(x_sph,y_sph,z_sph,'FaceColor','none','EdgeColor',lightGrey)
%         axis([-6371 6371 -6371 6371 -6371 6371])
%         view(142.5,30)
%         grid on
%         hold on
%         tetramesh(EL2NOD4(:,115)',GCOORD','FaceColor',lightWhi,'FaceAlpha',0.3)
%         scatter3(cart_mesh_corrected(1,1), ...
%                  cart_mesh_corrected(2,1), ...
%                  cart_mesh_corrected(3,1),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_1)
%         scatter3(cart_mesh_corrected(1,2), ...
%                  cart_mesh_corrected(2,2), ...
%                  cart_mesh_corrected(3,2),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_2)
%         scatter3(cart_mesh_corrected(1,3), ...
%                  cart_mesh_corrected(2,3), ...
%                  cart_mesh_corrected(3,3),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_3)
%         scatter3(cart_mesh_corrected(1,4), ...
%                  cart_mesh_corrected(2,4), ...
%                  cart_mesh_corrected(3,4),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_4)
%         
%         scatter3(cart_mesh_corrected(1,5), ...
%                  cart_mesh_corrected(2,5), ...
%                  cart_mesh_corrected(3,5),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_1)
%         scatter3(cart_mesh_corrected(1,6), ...
%                  cart_mesh_corrected(2,6), ...
%                  cart_mesh_corrected(3,6),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_2)
%         scatter3(cart_mesh_corrected(1,7), ...
%                  cart_mesh_corrected(2,7), ...
%                  cart_mesh_corrected(3,7),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_2)
%         scatter3(cart_mesh_corrected(1,8), ...
%                  cart_mesh_corrected(2,8), ...
%                  cart_mesh_corrected(3,8),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_3)
%         scatter3(cart_mesh_corrected(1,9), ...
%                  cart_mesh_corrected(2,9), ...
%                  cart_mesh_corrected(3,9),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_3)
%         scatter3(cart_mesh_corrected(1,10), ...
%                  cart_mesh_corrected(2,10), ...
%                  cart_mesh_corrected(3,10),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_4)
%         scatter3(cart_mesh_corrected(1,11), ...
%                  cart_mesh_corrected(2,11), ...
%                  cart_mesh_corrected(3,11),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_4)
%         scatter3(cart_mesh_corrected(1,12), ...
%                  cart_mesh_corrected(2,12), ...
%                  cart_mesh_corrected(3,12),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_1)
%         scatter3(cart_mesh_corrected(1,13), ...
%                  cart_mesh_corrected(2,13), ...
%                  cart_mesh_corrected(3,13),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_1)
%         scatter3(cart_mesh_corrected(1,14), ...
%                  cart_mesh_corrected(2,14), ...
%                  cart_mesh_corrected(3,14),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_3)
%         scatter3(cart_mesh_corrected(1,15), ...
%                  cart_mesh_corrected(2,15), ...
%                  cart_mesh_corrected(3,15),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_4)
%         scatter3(cart_mesh_corrected(1,16), ...
%                  cart_mesh_corrected(2,16), ...
%                  cart_mesh_corrected(3,16),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_2)
%         title('element mesh corrected)')
%     end
    
    if any(i_out_cone_cross_2pi)
        % Check local node numbering in edges ouside the cone and crossing 2pi
        edges_out_cone_cross_2pi       = find(i_out_cone_cross_2pi > 0);
        nodes_edges_out_cone_cross_2pi = [2*edges_out_cone_cross_2pi-1 2*edges_out_cone_cross_2pi]'; 
        nodes_edges_out_cone_cross_2pi = nodes_edges_out_cone_cross_2pi(:); % it needs to be a colum vector
        
        diffth_rot180 = abs(thn_check_rot180(nodes_edges_out_cone_cross_2pi,i) - thn_mesh_rot180(nodes_edges_out_cone_cross_2pi,i)) > 1e-12;
        diffph_rot180 = abs(phn_check_rot180(nodes_edges_out_cone_cross_2pi,i) - phn_mesh_rot180(nodes_edges_out_cone_cross_2pi,i)) > 1e-12;
        diffr_rot180  = abs(rn_check_rot180(nodes_edges_out_cone_cross_2pi,i)  - rn_mesh_rot180(nodes_edges_out_cone_cross_2pi,i))  > 1e-08;
        
        % Correct local node numbering in edges ouside the cone and crossing 2pi
        for inod = 1:2:size(nodes_edges_out_cone_cross_2pi,1)
            if (diffth_rot180(inod) || diffph_rot180(inod) || diffr_rot180(inod))
                EL2NOD20([4+nodes_edges_out_cone_cross_2pi(inod) 4+nodes_edges_out_cone_cross_2pi(inod+1)],i) = ...
                    EL2NOD20([4+nodes_edges_out_cone_cross_2pi(inod+1) 4+nodes_edges_out_cone_cross_2pi(inod)],i);
            end
        end
    end
end


%=============================================================================================================================================
% COMPUTE NODES AT 1/3 and 2/3 FOR EDGES THAT ARE CROSSING THE CONE BOUNDARY IN THE 90° ROTATED FRAME AROUND X AXIS 
%=============================================================================================================================================
% Select those edges crossing the cone boundary
EDGE2NOD_cross_cone = EDGE2NOD;
[~,LOCA]            = ismember(EDGE2NOD_out_cone,EDGE2NOD,'rows');
[~,LOCA2]           = ismember(fliplr(EDGE2NOD_out_cone),EDGE2NOD,'rows'); % in case some edges were flipped when using unique_keep_order
LOCA(LOCA == 0)     = LOCA2(LOCA2 ~= 0); clear LOCA2;
EDGE2NOD_cross_cone(LOCA,:) = []; % remove those ouside cone edges in the EDGE2NOD_cross_cone list

% COMPUTE COORDINATES FOR NEW NODES ON THE EDGES
% ----------------------------------------------
RR_X_90_CCW     = [ 1  0  0; 0  0 -1; 0  1  0 ];    % 90° counterclockwise around the X axis rotation matrix
GCOORD4_rot     = RR_X_90_CCW * GCOORD4;            % rotate coordinates
GCOORD4_SPH_rot = cartesian2spherical(GCOORD4_rot); % convert to spherical coordinates

% Coordinates of the nodes defining each bar's end points
nedge_cross_cone       = size(EDGE2NOD_cross_cone,1); % number of unique edges (i.e. number of new edge nodes that have to be generated)
thBarEnds_cross_cone   = reshape(GCOORD4_SPH_rot(1,EDGE2NOD_cross_cone'),2,nedge_cross_cone);
phBarEnds_cross_cone   = reshape(GCOORD4_SPH_rot(2,EDGE2NOD_cross_cone'),2,nedge_cross_cone);
rBarEnds_cross_cone    = reshape(GCOORD4_SPH_rot(3,EDGE2NOD_cross_cone'),2,nedge_cross_cone);

% Create two new nodes on each edge
LthBar_cross_cone      = thBarEnds_cross_cone(2,:) - thBarEnds_cross_cone(1,:);
thBarThirds_cross_cone = [thBarEnds_cross_cone(1,:)+(1/3)*LthBar_cross_cone
                          thBarEnds_cross_cone(1,:)+(2/3)*LthBar_cross_cone];
thBarThirds_cross_cone = thBarThirds_cross_cone(:)';

LphBar_els_cross_cone  = phBarEnds_cross_cone(2,:) - phBarEnds_cross_cone(1,:);
phBarThirds_cross_cone = [phBarEnds_cross_cone(1,:)+(1/3)*LphBar_els_cross_cone
                          phBarEnds_cross_cone(1,:)+(2/3)*LphBar_els_cross_cone];
phBarThirds_cross_cone = phBarThirds_cross_cone(:)';

LrBar_cross_cone       = rBarEnds_cross_cone(2,:) - rBarEnds_cross_cone(1,:);
rBarThirds_cross_cone  = [rBarEnds_cross_cone(1,:)+(1/3)*LrBar_cross_cone
                          rBarEnds_cross_cone(1,:)+(2/3)*LrBar_cross_cone];
rBarThirds_cross_cone  = rBarThirds_cross_cone(:)';

GCOORD_SPH_bars_cross_cone_rot = [thBarThirds_cross_cone; ...
                                  phBarThirds_cross_cone; ...
                                  rBarThirds_cross_cone];

GCOORD_bars_cross_cone_rot = spherical2cartesian(GCOORD_SPH_bars_cross_cone_rot); % convert to Cartesian coordinates
RR_X_90_CW                 = [ 1  0  0; 0  0  1; 0 -1  0 ];                       % 90° clockwise around the X axis rotation matrix
GCOORD_bars_cross_cone     = RR_X_90_CW * GCOORD_bars_cross_cone_rot;             % undo the rotation
GCOORD_SPH_bars_cross_cone = cartesian2spherical(GCOORD_bars_cross_cone);         % convert to spherical coodinates

%=========================================================================================================================
% MERGE NODES AT 1/3 and 2/3 COMPUTED IN BOTH SPH FRAMES:
%   - NODES AT 1/3 and 2/3 OUTSIDE THE CONE BOUNDARY (ORIGINAL FRAME)
%   - NODES AT 1/3 and 2/3 CROSSING THE CONE BOUNDARY (90° ROTATED FRAME)
%=========================================================================================================================
if ~isempty(EDGE2NOD_cross_cone)
    % Compute positions of EDGE2NOD_cross_cone in the EDGE2NOD list
    [~,LOCB]        = ismember(EDGE2NOD_cross_cone,EDGE2NOD,'rows');
    [~,LOCB2]       = ismember(fliplr(EDGE2NOD_cross_cone),EDGE2NOD,'rows'); % this is just in case some bars were flipped when using find_unique_bars
    LOCB(LOCB == 0) = LOCB2(LOCB2 ~= 0); %clear LOCB2;
    LOCD = [(LOCB*2-1)'; (LOCB*2)'];
    LOCD = LOCD(:)';
    % Substitute nodes of edges crossing the cone
    thBarThirds(LOCD) = GCOORD_SPH_bars_cross_cone(1,:);
    phBarThirds(LOCD) = GCOORD_SPH_bars_cross_cone(2,:);
    rBarThirds(LOCD)  = GCOORD_SPH_bars_cross_cone(3,:);
end

% storage for vertices (4) and nodes on the edges (12) in original frame
GCOORD16_SPH                   = zeros(3,nnod16);
GCOORD16_SPH(:,1:nnod4)        = GCOORD4_SPH;
GCOORD16_SPH(1,nnod4+1:nnod16) = thBarThirds;
GCOORD16_SPH(2,nnod4+1:nnod16) = phBarThirds;
GCOORD16_SPH(3,nnod4+1:nnod16) = rBarThirds;
GCOORD16                       = spherical2cartesian(GCOORD16_SPH); % convert to Cartesian coordinates

%=========================================================================================================================
% Now correct local node numbering for edge nodes (crossing the cone boundary) 
%=========================================================================================================================
% Coordinates for vertices (4) and nodes on the edges (12) in the 90° rotated around X axis frame
GCOORD_rot90     = RR_X_90_CCW * GCOORD16;            % rotate coordinates
GCOORD_SPH_rot90 = cartesian2spherical(GCOORD_rot90); % convert to spherical coordinates

% Compute (theta,phi,r) given by the mesh coordinates and (theta,phi,r) given by the shape functions in the 180° rotated around Z axis frame 
thn_vertx_rot90  = reshape(GCOORD_SPH_rot90(1,EL2NOD20(1:nvertx,:)),nvertx,[]);
thn_check_rot90  = N' * thn_vertx_rot90;
thn_mesh_rot90   = reshape(GCOORD_SPH_rot90(1,EL2NOD20(5:16,:)),12,[]);
phn_vertx_rot90  = reshape(GCOORD_SPH_rot90(2,EL2NOD20(1:nvertx,:)),nvertx,[]);
phn_check_rot90  = N' * phn_vertx_rot90;
phn_mesh_rot90   = reshape(GCOORD_SPH_rot90(2,EL2NOD20(5:16,:)),12,[]);
rn_vertx_rot90   = reshape(GCOORD_SPH_rot90(3,EL2NOD20(1:nvertx,:)),nvertx,[]);
rn_check_rot90   = N' * rn_vertx_rot90;
rn_mesh_rot90    = reshape(GCOORD_SPH_rot90(3,EL2NOD20(5:16,:)),12,[]);

% Compare coordinates for each element
for i=1:nel
    this_el           = EL2NOD4(:,i);
    edges_this_el     = reshape( this_el(edges',:),2,[] )'; % edges of this element
    i_cross_cone      = ismember(edges_this_el,EDGE2NOD_cross_cone,'rows');
    i_cross_cone_flip = ismember(fliplr(edges_this_el),EDGE2NOD_cross_cone,'rows');
    i_cross_cone      = (i_cross_cone + i_cross_cone_flip) > 0;
    
    if any(i_cross_cone)
        % Check local node numbering in edges ouside the cone and crossing 2pi
        edges_cross_cone       = find(i_cross_cone > 0);
        nodes_edges_cross_cone = [2*edges_cross_cone-1 2*edges_cross_cone]'; 
        nodes_edges_cross_cone = nodes_edges_cross_cone(:); % it needs to be a colum vector
        
        diffth_rot90 = abs(thn_check_rot90(nodes_edges_cross_cone,i) - thn_mesh_rot90(nodes_edges_cross_cone,i)) > 1e-12;
        diffph_rot90 = abs(phn_check_rot90(nodes_edges_cross_cone,i) - phn_mesh_rot90(nodes_edges_cross_cone,i)) > 1e-12;
        diffr_rot90  = abs(rn_check_rot90(nodes_edges_cross_cone,i)  - rn_mesh_rot90(nodes_edges_cross_cone,i))  > 1e-08;
        
        % Correct local node numbering in edges ouside the cone and crossing 2pi
        for inod = 1:2:size(nodes_edges_cross_cone,1)
            if (diffth_rot90(inod) || diffph_rot90(inod) || diffr_rot90(inod))
                EL2NOD20([4+nodes_edges_cross_cone(inod) 4+nodes_edges_cross_cone(inod+1)],i) = ...
                    EL2NOD20([4+nodes_edges_cross_cone(inod+1) 4+nodes_edges_cross_cone(inod)],i);
            end
        end
    end
end

%=========================================================================================================================
% COMPUTE ALL FACE-NODES IN THE ORIGINAL FRAME
%=========================================================================================================================
nnodel             = size(EL2NOD20,1); % nodes per element
nel                = size(EL2NOD20,2); % number of isoparametric elements
nnod0              = max(EL2NOD20(:));             % total number of nodes
% Define faces using 9 nodes (according to sf_dsf_tet.m) since for
% isoparametric elements some edges are curved and some other are straight
faces = [2 3 4  7  8  9 10 15 16 ... % face 1
         1 3 4 13 14  9 10 11 12 ... % face 2
         1 2 4  5  6 16 15 11 12 ... % face 3
         1 2 3  5  6  7  8 14 13];   % face 4

FACE2NOD_all_faces = reshape( EL2NOD20(faces',:),9,[] )';

% Find the faces that are shared by neighboring elements and return a unique list of faces.
[FACE2NOD_all_faces,~,ib] = unique_keep_order(FACE2NOD_all_faces);
nface                     = size(FACE2NOD_all_faces,1);

th_all_face_nodes = sum(reshape(GCOORD16_SPH(1,FACE2NOD_all_faces'),9,nface))./9;
ph_all_face_nodes = sum(reshape(GCOORD16_SPH(2,FACE2NOD_all_faces'),9,nface))./9;
r_all_face_nodes  = sum(reshape(GCOORD16_SPH(3,FACE2NOD_all_faces'),9,nface))./9;

% element to bar connectivity after doubles have been merged (removed)
EL2FACE = reshape(uint32(ib),4,nel);

% write face node connectivity for isoparametric elements
EL2NOD20(nnodel+1:nnodel+4,:) = nnod0 + EL2FACE;

%======================================================================================================================================
% COMPUTE FACE-NODES THAT ARE OUTSIDE THE CONE BOUNDARY AND CROSSING PHI = 2PI IN THE 180° ROTATED FRAME AROUND Z AXIS 
% (WE USE ONLY THE 3 VERTICES TO COMPUTE THE FACE-NODES SINCE THESE FACES HAVE STRAIGHT EDGES IN THE 180° ROTATED FRAME AROUND Z AXIS) 
%======================================================================================================================================
TH2FACES          = reshape(GCOORD4_SPH(1,FACE2NOD_all_faces(:,1:3)),[],3);
vertices_out_cone = sum(TH2FACES > SETTINGS.theta_cone*pi/180 & ...
                        TH2FACES < (180-SETTINGS.theta_cone)*pi/180,2); % vector for faces having x vertices outside the cone
                        % 1 --> one vertex is outside the cone boundary (this is not an isoparametric element)
                        % 2 --> two vertices are outside the cone boundary (no face outside the cone)
                        % 3 --> three vertices are outside the cone boundary (1 face outside the cone)
FACE2NOD_1_face_out_cone = FACE2NOD_all_faces(vertices_out_cone == 3,:);
PH2FACES          = reshape(GCOORD4_SPH(2,FACE2NOD_1_face_out_cone(:,1:3)),[],3);
ang_dist          = zeros(3,size(FACE2NOD_1_face_out_cone,1)); 
ang_dist(1,:)     = abs(PH2FACES(:,1)-PH2FACES(:,2)); % angular distance between vertex 1 and vertex 2
ang_dist(2,:)     = abs(PH2FACES(:,2)-PH2FACES(:,3)); % angular distance between vertex 2 and vertex 3
ang_dist(3,:)     = abs(PH2FACES(:,3)-PH2FACES(:,1)); % angular distance between vertex 3 and vertex 1
ang_dist_longer_than_pi = ang_dist > pi; % angular distance longer than pi 
                                         % (means that some nodes are in the 1st/5th quadrant and some others nodes are in the 4th/8th quadrant)
faces_cross_2pi                    = sum(ang_dist_longer_than_pi,1) > 0; % boolean vector for elements crossing phi = 2pi
FACE2NOD_1_face_out_cone_cross_2pi = FACE2NOD_1_face_out_cone(faces_cross_2pi,:);

RR_Z_180_CCW                    = [-1  0  0 ; 0 -1  0 ; 0  0  1 ];     % 180° counterclockwise around the Z axis rotation matrix
GCOORD_rot_180                  = RR_Z_180_CCW * GCOORD16;               % rotate coordinates
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
%=========================================================================================================================
% Remove those ouside cone faces in the FACE2NOD_faces_crossing_cone list taking into account the possible permutations
face_perms           = perms([3 2 1]);
[~,LOCB]             = ismember(FACE2NOD_1_face_out_cone(:,face_perms(1,:)),FACE2NOD_all_faces(:,1:3),'rows');
for i = 2:6
    [~,LOCB2]        = ismember(FACE2NOD_1_face_out_cone(:,face_perms(i,:)),FACE2NOD_all_faces(:,1:3),'rows');
    LOCB(LOCB2 ~= 0) = LOCB2(LOCB2 ~= 0);
end
FACE2NOD_faces_crossing_cone         = FACE2NOD_all_faces;
FACE2NOD_faces_crossing_cone(LOCB,:) = []; % remove those ouside cone faces in the FACE2NOD_faces_crossing_cone list

GCOORD_rot_90              = RR_X_90_CCW * GCOORD16;             % rotate coordinates
GCOORD_SPH_rot_90          = cartesian2spherical(GCOORD_rot_90); % convert to spherical coordinates
nface_crossing_cone        = size(FACE2NOD_faces_crossing_cone,1);
th_faces_crossing_cone_rot = sum(reshape(GCOORD_SPH_rot_90(1,FACE2NOD_faces_crossing_cone'),9,nface_crossing_cone))./9;
ph_faces_crossing_cone_rot = sum(reshape(GCOORD_SPH_rot_90(2,FACE2NOD_faces_crossing_cone'),9,nface_crossing_cone))./9;
r_faces_crossing_cone_rot  = sum(reshape(GCOORD_SPH_rot_90(3,FACE2NOD_faces_crossing_cone'),9,nface_crossing_cone))./9;

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
GCOORD20              = [GCOORD16 GCOORD_face_nodes];

%=========================================================================================================================
% OUTPUT DATA
%=========================================================================================================================
nnod10                  = size(GCOORD,2);
EL2NOD_cubic            = EL2NOD20;
EL2NOD_cubic(5:20,:)    = (nnod10-nnod4) + EL2NOD_cubic(5:20,:);

GCOORD_cubic            = GCOORD20;
GCOORD_cubic(:,1:nnod4) = [];

% %=========================================================================================================================
% % PLOT
% %=========================================================================================================================
% figure(741);clf
% GCOORD_plot = [GCOORD GCOORD_cubic];
% els_plot    = 115; %13; %75;
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
% tetramesh(EL2NOD_cubic(1:4,els_plot)',GCOORD','FaceColor',lightWhi,'FaceAlpha',0.3)
% color_local_vertex_1 = [1 0 0]; % red
% color_local_vertex_2 = [1 1 0]; % yellow
% color_local_vertex_3 = [0 1 0]; % green
% color_local_vertex_4 = [0 0 1]; % blue
% scatter3(GCOORD(1,EL2NOD_cubic(1,els_plot)), ...
%          GCOORD(2,EL2NOD_cubic(1,els_plot)), ...
%          GCOORD(3,EL2NOD_cubic(1,els_plot)),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_1)
% scatter3(GCOORD(1,EL2NOD_cubic(2,els_plot)), ...
%          GCOORD(2,EL2NOD_cubic(2,els_plot)), ...
%          GCOORD(3,EL2NOD_cubic(2,els_plot)),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_2)
% scatter3(GCOORD(1,EL2NOD_cubic(3,els_plot)), ...
%          GCOORD(2,EL2NOD_cubic(3,els_plot)), ...
%          GCOORD(3,EL2NOD_cubic(3,els_plot)),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_3)
% scatter3(GCOORD(1,EL2NOD_cubic(4,els_plot)), ...
%          GCOORD(2,EL2NOD_cubic(4,els_plot)), ...
%          GCOORD(3,EL2NOD_cubic(4,els_plot)),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_4)
% 
% scatter3(GCOORD_plot(1,EL2NOD_cubic(5,els_plot)), ...
%          GCOORD_plot(2,EL2NOD_cubic(5,els_plot)), ...
%          GCOORD_plot(3,EL2NOD_cubic(5,els_plot)),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_1)
% scatter3(GCOORD_plot(1,EL2NOD_cubic(6,els_plot)), ...
%          GCOORD_plot(2,EL2NOD_cubic(6,els_plot)), ...
%          GCOORD_plot(3,EL2NOD_cubic(6,els_plot)),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_2)
% scatter3(GCOORD_plot(1,EL2NOD_cubic(7,els_plot)), ...
%          GCOORD_plot(2,EL2NOD_cubic(7,els_plot)), ...
%          GCOORD_plot(3,EL2NOD_cubic(7,els_plot)),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_2)
% scatter3(GCOORD_plot(1,EL2NOD_cubic(8,els_plot)), ...
%          GCOORD_plot(2,EL2NOD_cubic(8,els_plot)), ...
%          GCOORD_plot(3,EL2NOD_cubic(8,els_plot)),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_3)
% scatter3(GCOORD_plot(1,EL2NOD_cubic(9,els_plot)), ...
%          GCOORD_plot(2,EL2NOD_cubic(9,els_plot)), ...
%          GCOORD_plot(3,EL2NOD_cubic(9,els_plot)),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_3)
% scatter3(GCOORD_plot(1,EL2NOD_cubic(10,els_plot)), ...
%          GCOORD_plot(2,EL2NOD_cubic(10,els_plot)), ...
%          GCOORD_plot(3,EL2NOD_cubic(10,els_plot)),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_4)
% scatter3(GCOORD_plot(1,EL2NOD_cubic(11,els_plot)), ...
%          GCOORD_plot(2,EL2NOD_cubic(11,els_plot)), ...
%          GCOORD_plot(3,EL2NOD_cubic(11,els_plot)),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_4)
% scatter3(GCOORD_plot(1,EL2NOD_cubic(12,els_plot)), ...
%          GCOORD_plot(2,EL2NOD_cubic(12,els_plot)), ...
%          GCOORD_plot(3,EL2NOD_cubic(12,els_plot)),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_1)
% scatter3(GCOORD_plot(1,EL2NOD_cubic(13,els_plot)), ...
%          GCOORD_plot(2,EL2NOD_cubic(13,els_plot)), ...
%          GCOORD_plot(3,EL2NOD_cubic(13,els_plot)),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_1)
% scatter3(GCOORD_plot(1,EL2NOD_cubic(14,els_plot)), ...
%          GCOORD_plot(2,EL2NOD_cubic(14,els_plot)), ...
%          GCOORD_plot(3,EL2NOD_cubic(14,els_plot)),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_3)
% scatter3(GCOORD_plot(1,EL2NOD_cubic(15,els_plot)), ...
%          GCOORD_plot(2,EL2NOD_cubic(15,els_plot)), ...
%          GCOORD_plot(3,EL2NOD_cubic(15,els_plot)),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_4)
% scatter3(GCOORD_plot(1,EL2NOD_cubic(16,els_plot)), ...
%          GCOORD_plot(2,EL2NOD_cubic(16,els_plot)), ...
%          GCOORD_plot(3,EL2NOD_cubic(16,els_plot)),'MarkerEdgeColor','k','MarkerFaceColor',color_local_vertex_2)
%      
% scatter3(GCOORD_plot(1,EL2NOD_cubic(17,els_plot)), ...
%          GCOORD_plot(2,EL2NOD_cubic(17,els_plot)), ...
%          GCOORD_plot(3,EL2NOD_cubic(17,els_plot)),'MarkerEdgeColor',color_local_vertex_1,'MarkerFaceColor',color_local_vertex_1)
% scatter3(GCOORD_plot(1,EL2NOD_cubic(18,els_plot)), ...
%          GCOORD_plot(2,EL2NOD_cubic(18,els_plot)), ...
%          GCOORD_plot(3,EL2NOD_cubic(18,els_plot)),'MarkerEdgeColor',color_local_vertex_2,'MarkerFaceColor',color_local_vertex_2)
% scatter3(GCOORD_plot(1,EL2NOD_cubic(19,els_plot)), ...
%          GCOORD_plot(2,EL2NOD_cubic(19,els_plot)), ...
%          GCOORD_plot(3,EL2NOD_cubic(19,els_plot)),'MarkerEdgeColor',color_local_vertex_3,'MarkerFaceColor',color_local_vertex_3)
% scatter3(GCOORD_plot(1,EL2NOD_cubic(20,els_plot)), ...
%          GCOORD_plot(2,EL2NOD_cubic(20,els_plot)), ...
%          GCOORD_plot(3,EL2NOD_cubic(20,els_plot)),'MarkerEdgeColor',color_local_vertex_4,'MarkerFaceColor',color_local_vertex_4)

end % END OF FUNCTION cubic_els_in_cone_iso

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