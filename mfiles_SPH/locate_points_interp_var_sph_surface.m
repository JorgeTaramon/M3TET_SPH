function [U_SPH_i,lc] = locate_points_interp_var_sph_surface(GCOORD_SPH,EL2NOD_face,els,gTH_PT,MESH,els_face,U_SPH)
% Usage: lc = local_coords_3d_sph_surface(GCOORD_SPH,EL2NOD,els,gTH_PT,MESH,els_face)
%
% Purpose:
%   Returns local coordinates of points "gTH_PT" in elements "els".
%
% Input:
%   GCOORD_SPH : [matrix]    : spherical coodinates (theta,phi,r) 
%                              (in rad, rad and km) of the mesh
%   EL2NOD     : [matrix]    : connectivity matrix for GCOORD_SPH
%   els        : [rowvector] : element in which each point is located
%   gTH_PT     : [matrix]    : coordinates of points to be located
%   MESH       : [structure] : FE mesh parameters
%
% Output:
%   lc         : [matrix]    : local coordinates of points in each element
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH March 2011
% JH March 2013
% JMT Oct 2016 : works for spherical coordinates

%===============================================================================================================================================================================================
% LOCAL COORDINATES FOR THOSE POINTS (gTH_PT) INSIDE ELEMENTS THAT ARE OUTSIDE THE CONE AND NOT CROSSING PHI = 2PI (STRAIGHT EDGES IN THE ORIGINAL SPHERICAL FRAME) 
%===============================================================================================================================================================================================
els_globlal = els_face(els);
[~,els_out_cone_no_cross_2pi]     = ismember(els_globlal(ismember(els_globlal,MESH.els_out_cone_no_cross_2pi)),els_face); % elements outside the cone and not crossing 2pi from elements in which each point is located (els)
gTH_PT_els_out_cone_no_cross_2pi  = gTH_PT(:,ismember(els,els_out_cone_no_cross_2pi)); % coordinates of points (inside elements outside the cone and not crossing 2pi) to be located
% correct local node numbering
EL2NOD_face                       = correct_local_node_numbering(GCOORD_SPH,EL2NOD_face,els_out_cone_no_cross_2pi);
% compute local coordinates
lc_els_out_cone_no_cross_2pi      = compute_lc(GCOORD_SPH(1:2,:),EL2NOD_face,els_out_cone_no_cross_2pi,gTH_PT_els_out_cone_no_cross_2pi); % *SUBFUNCTION*
% DShape function values at local coordinates
N   = sf_dsf_tri36712(lc_els_out_cone_no_cross_2pi,3,'matrix');
% Interpolate temperature to points using shape functions
Uth                               = U_SPH(1,:);
Uph                               = U_SPH(2,:);
Uth_i_els_out_cone_no_cross_2pi   = sum(N.*Uth(EL2NOD_face(1:3,els_out_cone_no_cross_2pi)));
Uph_i_els_out_cone_no_cross_2pi   = sum(N.*Uph(EL2NOD_face(1:3,els_out_cone_no_cross_2pi)));
U_SPH_i_els_out_cone_no_cross_2pi = [Uth_i_els_out_cone_no_cross_2pi; Uph_i_els_out_cone_no_cross_2pi; zeros(1,size(els_out_cone_no_cross_2pi,2))];

figure(990);clf
plot(gTH_PT_els_out_cone_no_cross_2pi(1,:),gTH_PT_els_out_cone_no_cross_2pi(2,:),'xr')
axis equal
hold on
trimesh(EL2NOD_face(:,els_out_cone_no_cross_2pi)',GCOORD_SPH(1,:)',GCOORD_SPH(2,:)','Color',[0 0 0])
xlabel('\theta (rad)')
ylabel('\phi (rad)')

%===============================================================================================================================================================================================
% LOCAL COORDINATES FOR THOSE POINTS (gTH_PT) INSIDE ELEMENTS THAT ARE OUTSIDE THE CONE AND CROSSING PHI = 2PI (STRAIGHT EDGES IN THE SPHERICAL ROTATED FRAME (180° AROUND Z AXIS)) 
%===============================================================================================================================================================================================
[~,els_out_cone_cross_2pi]        = ismember(els_globlal(ismember(els_globlal,MESH.els_out_cone_cross_2pi)),els_face); % elements outside the cone and crossing 2pi from elements in which each point is located (els)
gTH_PT_els_out_cone_cross_2pi     = gTH_PT(:,ismember(els,els_out_cone_cross_2pi)); % coordinates of points (inside elements outside the cone and crossing 2pi) to be located
gX_PT_els_out_cone_cross_2pi      = spherical2cartesian(gTH_PT_els_out_cone_cross_2pi);
RR_Z_180_CCW                      = [-1  0  0 ; ...
                                      0 -1  0 ; ...
                                      0  0  1 ];   % rotation matrix
gX_PT_els_out_cone_cross_2pi_rot  = RR_Z_180_CCW * gX_PT_els_out_cone_cross_2pi;           % rotate those points
gTH_PT_els_out_cone_cross_2pi_rot = cartesian2spherical(gX_PT_els_out_cone_cross_2pi_rot); % convert back to spherical coordinates
GCOORD_rot_180_Z                  = RR_Z_180_CCW * MESH.GCOORD;                            % rotate all mesh nodes in Cartesian coordinates
GCOORD_SPH_rot_180_Z              = cartesian2spherical(GCOORD_rot_180_Z);                 % convert back to spherical coordinates
% correct local node numbering
EL2NOD_face                       = correct_local_node_numbering(GCOORD_SPH_rot_180_Z,EL2NOD_face,els_out_cone_cross_2pi);
% compute local coordinates
lc_els_out_cone_cross_2pi         = compute_lc(GCOORD_SPH_rot_180_Z(1:2,:),EL2NOD_face,els_out_cone_cross_2pi,gTH_PT_els_out_cone_cross_2pi_rot); % *SUBFUNCTION*
% DShape function values at local coordinates
N                                 = sf_dsf_tri36712(lc_els_out_cone_cross_2pi,3,'matrix');
% Interpolate temperature to points using shape functions
TT1                               = make_rotation_matrix(MESH.GCOORD);
U                                 = TT1 * U_SPH(:);
U_rot                             = RR_Z_180_CCW * reshape(U,3,[]);
TT2                               = make_rotation_matrix(GCOORD_rot_180_Z);
U_SPH_rot                         = TT2' * U_rot(:);
U_SPH_rot                         = reshape(U_SPH_rot,3,[]);
Uth                               = U_SPH_rot(1,:);
Uph                               = U_SPH_rot(2,:);
Uth_i_els_out_cone_cross_2pi      = sum(N.*Uth(EL2NOD_face(1:3,els_out_cone_cross_2pi)));
Uph_i_els_out_cone_cross_2pi      = sum(N.*Uph(EL2NOD_face(1:3,els_out_cone_cross_2pi)));
U_SPH_rot_i                       = [Uth_i_els_out_cone_cross_2pi; Uph_i_els_out_cone_cross_2pi; zeros(1,size(els_out_cone_cross_2pi,2))];
TT3                               = make_rotation_matrix(gX_PT_els_out_cone_cross_2pi_rot);
U_rot_i                           = TT3 * U_SPH_rot_i(:);
U_i                               = RR_Z_180_CCW' * reshape(U_rot_i,3,[]);
TT4                               = make_rotation_matrix(gX_PT_els_out_cone_cross_2pi);
U_SPH_i                           = TT4' * U_i(:);
U_SPH_i_els_out_cone_cross_2pi    = reshape(U_SPH_i,3,[]);

figure(991);clf
plot(gTH_PT_els_out_cone_cross_2pi_rot(1,:),gTH_PT_els_out_cone_cross_2pi_rot(2,:),'xr')
axis equal
hold on
trimesh(EL2NOD_face(:,els_out_cone_cross_2pi)',GCOORD_SPH_rot_180_Z(1,:)',GCOORD_SPH_rot_180_Z(2,:)','Color',[0 0 0])
xlabel('\theta (rad)')
ylabel('\phi (rad)')

%===============================================================================================================================================================================================
% LOCAL COORDINATES FOR THOSE POINTS (gTH_PT) INSIDE ELEMENTS THAT ARE WITHIN THE CONE AND NOT ISOPARAMETRIC (STRAIGHT EDGES IN THE SPHERICAL ROTATED FRAME (90° AROUND X AXIS)) 
%===============================================================================================================================================================================================
[~,els_in_cone_no_iso]            = ismember(els_globlal(ismember(els_globlal,MESH.els_in_cone_no_iso)),els_face); % elements inside the cone and not isoparametric from elements in which each point is located (els)
gTH_PT_els_in_cone_no_iso         = gTH_PT(:,ismember(els,els_in_cone_no_iso)); % coordinates of points (inside elements within the cone and not isoparametric) to be located
gX_PT_els_in_cone_no_iso          = spherical2cartesian(gTH_PT_els_in_cone_no_iso);
RR_X_90_CCW                       = [ 1  0  0 ; ...
                                      0  0 -1 ; ...
                                      0  1  0 ];   % rotation matrix
gX_PT_els_in_cone_no_iso_rot      = RR_X_90_CCW * gX_PT_els_in_cone_no_iso;            % rotate those points
gTH_PT_els_in_cone_no_iso_rot     = cartesian2spherical(gX_PT_els_in_cone_no_iso_rot); % convert back to spherical coordinates
GCOORD_rot_90_X                   = RR_X_90_CCW * MESH.GCOORD;                         % rotate all mesh nodes in Cartesian coordinates
GCOORD_SPH_rot_90_X               = cartesian2spherical(GCOORD_rot_90_X);              % convert back to spherical coordinates
% correct local node numbering
EL2NOD_face                       = correct_local_node_numbering(GCOORD_SPH_rot_90_X,EL2NOD_face,els_in_cone_no_iso);
% compute local coordinates
lc_els_in_cone_no_iso             = compute_lc(GCOORD_SPH_rot_90_X(1:2,:),EL2NOD_face,els_in_cone_no_iso,gTH_PT_els_in_cone_no_iso_rot); % *SUBFUNCTION*
% DShape function values at local coordinates
N   = sf_dsf_tri36712(lc_els_in_cone_no_iso,3,'matrix');
% Interpolate temperature to points using shape functions
TT1                               = make_rotation_matrix(MESH.GCOORD);
U                                 = TT1 * U_SPH(:);
U_rot                             = RR_X_90_CCW * reshape(U,3,[]);
TT2                               = make_rotation_matrix(GCOORD_rot_90_X);
U_SPH_rot                         = TT2' * U_rot(:);
U_SPH_rot                         = reshape(U_SPH_rot,3,[]);
Uth                               = U_SPH_rot(1,:);
Uph                               = U_SPH_rot(2,:);
Uth_i_els_in_cone_no_iso          = sum(N.*Uth(EL2NOD_face(1:3,els_in_cone_no_iso)));
Uph_i_els_in_cone_no_iso          = sum(N.*Uph(EL2NOD_face(1:3,els_in_cone_no_iso)));
U_SPH_rot_i                       = [Uth_i_els_in_cone_no_iso; Uph_i_els_in_cone_no_iso; zeros(1,size(els_in_cone_no_iso,2))];
TT3                               = make_rotation_matrix(gX_PT_els_in_cone_no_iso_rot);
U_rot_i                           = TT3 * U_SPH_rot_i(:);
U_i                               = RR_X_90_CCW' * reshape(U_rot_i,3,[]);
TT4                               = make_rotation_matrix(gX_PT_els_in_cone_no_iso);
U_SPH_i                           = TT4' * U_i(:);
U_SPH_i_els_in_cone_no_iso        = reshape(U_SPH_i,3,[]);

figure(992);clf
plot(gTH_PT_els_in_cone_no_iso_rot(1,:),gTH_PT_els_in_cone_no_iso_rot(2,:),'xr')
axis equal
hold on
trimesh(EL2NOD_face(:,els_in_cone_no_iso)',GCOORD_SPH_rot_90_X(1,:)',GCOORD_SPH_rot_90_X(2,:)','Color',[0 0 0])
xlabel('\theta (rad)')
ylabel('\phi (rad)')

%===============================================================================================================================================================================================
% LOCAL COORDINATES FOR THOSE POINTS (gTH_PT) INSIDE ELEMENTS THAT ARE WITHIN THE CONE AND ISOPARAMETRIC (CURVED EDGES IN THE SPHERICAL ROTATED FRAME (90° AROUND X AXIS)) 
%===============================================================================================================================================================================================
[~,els_in_cone_iso]               = ismember(els_globlal(ismember(els_globlal,MESH.els_in_cone_iso)),els_face); % elements inside the cone and isoparametric from elements in which each point is located (els)
gTH_PT_els_in_cone_iso            = gTH_PT(:,ismember(els,els_in_cone_iso)); % coordinates of points (inside elements within the cone and isoparametric) to be located
gX_PT_els_in_cone_iso             = spherical2cartesian(gTH_PT_els_in_cone_iso);
RR_X_90_CCW                       = [ 1  0  0 ; ...
                                      0  0 -1 ; ...
                                      0  1  0 ];   % rotation matrix
gX_PT_els_in_cone_iso_rot         = RR_X_90_CCW * gX_PT_els_in_cone_iso;            % rotate those points
gTH_PT_els_in_cone_iso_rot        = cartesian2spherical(gX_PT_els_in_cone_iso_rot); % convert back to spherical coordinates
% correct local node numbering
EL2NOD_face                       = correct_local_node_numbering(GCOORD_SPH_rot_90_X,EL2NOD_face,els_in_cone_iso);
% compute local coordinates
lc_els_in_cone_iso                = compute_lc(GCOORD_SPH_rot_90_X(1:2,:),EL2NOD_face,els_in_cone_iso,gTH_PT_els_in_cone_iso_rot); % *SUBFUNCTION*
% lc_els_in_cone_iso                = compute_lc_surface_curved(GCOORD_SPH_rot_90_X(1:2,:),EL2NOD_face,els_in_cone_iso,gTH_PT_els_in_cone_iso_rot);
% DShape function values at local coordinates
N   = sf_dsf_tri36712(lc_els_in_cone_iso,3,'matrix');
% Interpolate temperature to points using shape functions
TT1                               = make_rotation_matrix(MESH.GCOORD);
U                                 = TT1 * U_SPH(:);
U_rot                             = RR_X_90_CCW * reshape(U,3,[]);
TT2                               = make_rotation_matrix(GCOORD_rot_90_X);
U_SPH_rot                         = TT2' * U_rot(:);
U_SPH_rot                         = reshape(U_SPH_rot,3,[]);
Uth                               = U_SPH_rot(1,:);
Uph                               = U_SPH_rot(2,:);
Uth_i_els_in_cone_iso             = sum(N.*Uth(EL2NOD_face(1:3,els_in_cone_iso)));
Uph_i_els_in_cone_iso             = sum(N.*Uph(EL2NOD_face(1:3,els_in_cone_iso)));
U_SPH_rot_i                       = [Uth_i_els_in_cone_iso; Uph_i_els_in_cone_iso; zeros(1,size(els_in_cone_iso,2))];
TT3                               = make_rotation_matrix(gX_PT_els_in_cone_iso_rot);
U_rot_i                           = TT3 * U_SPH_rot_i(:);
U_i                               = RR_X_90_CCW' * reshape(U_rot_i,3,[]);
TT4                               = make_rotation_matrix(gX_PT_els_in_cone_iso);
U_SPH_i                           = TT4' * U_i(:);
U_SPH_i_els_in_cone_iso           = reshape(U_SPH_i,3,[]);

figure(993);clf
plot(gTH_PT_els_in_cone_iso_rot(1,:),gTH_PT_els_in_cone_iso_rot(2,:),'xr')
axis equal
hold on
trimesh(EL2NOD_face(:,els_in_cone_iso)',GCOORD_SPH_rot_90_X(1,:)',GCOORD_SPH_rot_90_X(2,:)','Color',[0 0 0])
plot(GCOORD_SPH_rot_90_X(1,EL2NOD_face(:,els_in_cone_iso)),GCOORD_SPH_rot_90_X(2,EL2NOD_face(:,els_in_cone_iso)),'+k')
xlabel('\theta (rad)')
ylabel('\phi (rad)')

%===============================================================================================================================================================================================
% LOCAL COORDINATES FOR POINTS (gTH_PT)
%===============================================================================================================================================================================================
U_SPH_i                                                         = zeros(3,size(gTH_PT,2));
U_SPH_i(:,ismember(els_globlal,MESH.els_out_cone_no_cross_2pi)) = U_SPH_i_els_out_cone_no_cross_2pi;
U_SPH_i(:,ismember(els_globlal,MESH.els_out_cone_cross_2pi))    = U_SPH_i_els_out_cone_cross_2pi;
U_SPH_i(:,ismember(els_globlal,MESH.els_in_cone_no_iso))        = U_SPH_i_els_in_cone_no_iso;
U_SPH_i(:,ismember(els_globlal,MESH.els_in_cone_iso))           = U_SPH_i_els_in_cone_iso;

%===============================================================================================================================================================================================
% LOCAL COORDINATES FOR POINTS (gTH_PT)
%===============================================================================================================================================================================================
lc                                                         = zeros(2,size(gTH_PT,2));
lc(:,ismember(els_globlal,MESH.els_out_cone_no_cross_2pi)) = lc_els_out_cone_no_cross_2pi;
lc(:,ismember(els_globlal,MESH.els_out_cone_cross_2pi))    = lc_els_out_cone_cross_2pi;
lc(:,ismember(els_globlal,MESH.els_in_cone_no_iso))        = lc_els_in_cone_no_iso;
lc(:,ismember(els_globlal,MESH.els_in_cone_iso))           = lc_els_in_cone_iso;

end % END OF FUNCTION local_coords_3d_sph

% ##############################################################################################################################################################################################
%                               SUBFUNCTIONS
% ##############################################################################################################################################################################################

function EL2NOD_face = correct_local_node_numbering(GCOORD_SPH,EL2NOD_face,els)

nvertx = 3;
% Local node numbering:
%
%        3
%        | \
% s-axis 6   5
%        |     \
%        1 - 4 - 2
%          r-axis
%
% local coordinates of all vertices (1:3)
r = [0 1 0];
s = [0 0 1];
N = sf_dsf_tri36712([r;s],nvertx,'matrix');
thn_vertx = reshape(GCOORD_SPH(1,EL2NOD_face(1:nvertx,els)),nvertx,[]);
thn_check = N' * thn_vertx;
thn_mesh  = reshape(GCOORD_SPH(1,EL2NOD_face(1:nvertx,els)),nvertx,[]);
if ~isequal(thn_check,thn_mesh)
    error('wrong local numbering for vertices')
end
% local coordinates of all midside nodes (4:6)
r = [0.5 0.5 0  ];
s = [0   0.5 0.5];
N = sf_dsf_tri36712([r;s],nvertx,'matrix');
thn_vertx = reshape(GCOORD_SPH(1,EL2NOD_face(1:nvertx,els)),nvertx,[]);
thn_check = N' * thn_vertx;
thn_mesh  = reshape(GCOORD_SPH(1,EL2NOD_face(4:6,els)),3,[]);
diffth    = abs(thn_check-thn_mesh)>1e-12;
phn_vertx = reshape(GCOORD_SPH(2,EL2NOD_face(1:nvertx,els)),nvertx,[]);
phn_check = N' * phn_vertx;
phn_mesh  = reshape(GCOORD_SPH(2,EL2NOD_face(4:6,els)),3,[]);
diffph    = abs(phn_check-phn_mesh)>1e-12;

change_2_nodes = find(sum(diffth) == 2 | sum(diffph) == 2);
nod2flip = repmat((4:6)',1,size(change_2_nodes,2));
nod2flip = reshape(nod2flip(diffth(:,change_2_nodes)),2,[]);
for i=1:size(change_2_nodes,2)
    EL2NOD_face(nod2flip(:,i),els(change_2_nodes(i))) = EL2NOD_face(flipud(nod2flip(:,i)),els(change_2_nodes(i)));
end

change_3_nodes = find(sum(diffth) == 3 | sum(diffph) == 3);
for j=1:size(change_3_nodes,2)
    [I,nod2flip] = ismember(thn_check(:,change_3_nodes(j))',thn_mesh(:,change_3_nodes(j))');
    if sum(I)~= 3
        perm = [1 2 3; 2 3 1; 3 1 2; 1 3 2; 2 1 3; 3 2 1];
        a = thn_check(:,change_3_nodes(j))';
        b = thn_mesh(:,change_3_nodes(j))';
        c = zeros(6,1);
        for k=1:6
            c(k) = sum(abs(a-b(perm(k,:))));
        end
        nod2flip = perm(c == min(c),:);
    end
    nod2flip = 3 + nod2flip(1,:)';
    EL2NOD_face(4:6,els(change_3_nodes(j))) = EL2NOD_face(nod2flip,els(change_3_nodes(j)));
end

end % END OF SUBFUNCTION correct_local_node_numbering

% ##############################################################################################################################################################################################

function lc = compute_lc(GCOORD_SPH,EL2NOD_face,els,gTH_PT)
% Usage: lc = compute_lc(GCOORD_SPH,EL2NOD_face,els,gTH_PT)
%
% Purpose: Returns local coordinates of points "gTH_PT" in elements "els".
%
% Input:
%   GCOORD_SPH  : [matrix]    : spherical coodinates (theta,phi,r) 
%                               (in rad, rad and km) of the mesh 
%   EL2NOD_face : [matrix]    : connectivity matrix for faces on the
%                               surface
%   els         : [rowvector] : element in which each point is located
%   gTH_PT      : [matrix]    : coordinates of points to be located
%
% Output:
%   lc          : [matrix]    : local coordinates of points in each element
%
% Part of M2TRI - 2D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% Definition of local coordinates
% Make sure that shape (interpolation) functions are defined accordingly!
%
%   / \    3
%    |     | \
% lc(2,:)  |   \
%    |     |     \
%    |     1 - - - 2     
%        -- lc(1,:) -->
%
% JH Jan 2011
% JH Jan 2013 : handles NaN in els
% JH Nov 2014 : verified; added comment
% JMT May 2016 : works for polar coordinates

ind = ~isnan(els) & els>0;
els = els(ind);
ns  = length(els);
th  = reshape(GCOORD_SPH(1,EL2NOD_face(1:3,els)),3,ns);
ph  = reshape(GCOORD_SPH(2,EL2NOD_face(1:3,els)),3,ns);

thp = gTH_PT(1,ind);
php = gTH_PT(2,ind);

lc        = nan(2,length(ind));
denom     = -th(1,:).*ph(3,:)+th(1,:).*ph(2,:)-th(2,:).*ph(1,:)+th(2,:).*ph(3,:)+th(3,:).*ph(1,:)-th(3,:).*ph(2,:);
lc(1,ind) = -(-th(1,:).*php+th(1,:).*ph(3,:)-th(3,:).*ph(1,:)+thp.*ph(1,:)+th(3,:).*php-thp.*ph(3,:)) ./ denom;
lc(2,ind) =  (th(1,:).*ph(2,:)-th(1,:).*php-th(2,:).*ph(1,:)+th(2,:).*php+thp.*ph(1,:)-thp.*ph(2,:)) ./ denom;

end % % END OF FUNCTION compute_lc