function lc = local_coords_3d_sph(GCOORD_SPH,EL2NOD,els,gTH_PT,MESH)
% Usage: lc = local_coords_3d_sph(GCOORD_SPH,EL2NOD,els,gTH_PT,MESH)
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
els_out_cone_no_cross_2pi         = els(ismember(els,MESH.els_out_cone_no_cross_2pi)); % elements outside the cone and not crossing 2pi from elements in which each point is located (els)
gTH_PT_els_out_cone_no_cross_2pi  = gTH_PT(:,ismember(els,els_out_cone_no_cross_2pi)); % coordinates of points (inside elements outside the cone and not crossing 2pi) to be located
lc_els_out_cone_no_cross_2pi      = compute_lc(GCOORD_SPH,EL2NOD,els_out_cone_no_cross_2pi,gTH_PT_els_out_cone_no_cross_2pi); % *SUBFUNCTION*

%===============================================================================================================================================================================================
% LOCAL COORDINATES FOR THOSE POINTS (gTH_PT) INSIDE ELEMENTS THAT ARE OUTSIDE THE CONE AND CROSSING PHI = 2PI (STRAIGHT EDGES IN THE SPHERICAL ROTATED FRAME (180° AROUND Z AXIS)) 
%===============================================================================================================================================================================================
els_out_cone_cross_2pi            = els(ismember(els,MESH.els_out_cone_cross_2pi)); % elements outside the cone and crossing 2pi from elements in which each point is located (els)
gTH_PT_els_out_cone_cross_2pi     = gTH_PT(:,ismember(els,els_out_cone_cross_2pi)); % coordinates of points (inside elements outside the cone and crossing 2pi) to be located
gX_PT_els_out_cone_cross_2pi      = spherical2cartesian(gTH_PT_els_out_cone_cross_2pi);
RR_Z_180_CCW                      = [-1  0  0 ; ...
                                      0 -1  0 ; ...
                                      0  0  1 ];   % rotation matrix
gX_PT_els_out_cone_cross_2pi_rot  = RR_Z_180_CCW * gX_PT_els_out_cone_cross_2pi;           % rotate those points
gTH_PT_els_out_cone_cross_2pi_rot = cartesian2spherical(gX_PT_els_out_cone_cross_2pi_rot); % convert back to spherical coordinates
GCOORD_rot_180_Z                  = RR_Z_180_CCW * MESH.GCOORD;                            % rotate all mesh nodes in Cartesian coordinates
GCOORD_SPH_rot_180_Z              = cartesian2spherical(GCOORD_rot_180_Z);                 % convert back to spherical coordinates
lc_els_out_cone_cross_2pi         = compute_lc(GCOORD_SPH_rot_180_Z,EL2NOD,els_out_cone_cross_2pi,gTH_PT_els_out_cone_cross_2pi_rot); % *SUBFUNCTION*

%===============================================================================================================================================================================================
% LOCAL COORDINATES FOR THOSE POINTS (gTH_PT) INSIDE ELEMENTS THAT ARE WITHIN THE CONE AND NOT ISOPARAMETRIC (STRAIGHT EDGES IN THE SPHERICAL ROTATED FRAME (90° AROUND X AXIS)) 
%===============================================================================================================================================================================================
els_in_cone_no_iso                = els(ismember(els,MESH.els_in_cone_no_iso)); % elements inside the cone and not isoparametric from elements in which each point is located (els)
gTH_PT_els_in_cone_no_iso         = gTH_PT(:,ismember(els,els_in_cone_no_iso)); % coordinates of points (inside elements within the cone and not isoparametric) to be located
gX_PT_els_in_cone_no_iso          = spherical2cartesian(gTH_PT_els_in_cone_no_iso);
RR_X_90_CCW                       = [ 1  0  0 ; ...
                                      0  0 -1 ; ...
                                      0  1  0 ];   % rotation matrix
gX_PT_els_in_cone_no_iso_rot      = RR_X_90_CCW * gX_PT_els_in_cone_no_iso;            % rotate those points
gTH_PT_els_in_cone_no_iso_rot     = cartesian2spherical(gX_PT_els_in_cone_no_iso_rot); % convert back to spherical coordinates
GCOORD_rot_90_X                   = RR_X_90_CCW * MESH.GCOORD;                         % rotate all mesh nodes in Cartesian coordinates
GCOORD_SPH_rot_90_X               = cartesian2spherical(GCOORD_rot_90_X);              % convert back to spherical coordinates
lc_els_in_cone_no_iso             = compute_lc(GCOORD_SPH_rot_90_X,EL2NOD,els_in_cone_no_iso,gTH_PT_els_in_cone_no_iso_rot); % *SUBFUNCTION*

%===============================================================================================================================================================================================
% LOCAL COORDINATES FOR THOSE POINTS (gTH_PT) INSIDE ELEMENTS THAT ARE WITHIN THE CONE AND ISOPARAMETRIC (CURVED EDGES IN THE SPHERICAL ROTATED FRAME (90° AROUND X AXIS)) 
%===============================================================================================================================================================================================
els_in_cone_iso                   = els(ismember(els,MESH.els_in_cone_iso)); % elements inside the cone and isoparametric from elements in which each point is located (els)
gTH_PT_els_in_cone_iso            = gTH_PT(:,ismember(els,els_in_cone_iso)); % coordinates of points (inside elements within the cone and isoparametric) to be located
gX_PT_els_in_cone_iso             = spherical2cartesian(gTH_PT_els_in_cone_iso);
RR_X_90_CCW                       = [ 1  0  0 ; ...
                                      0  0 -1 ; ...
                                      0  1  0 ];   % rotation matrix
gX_PT_els_in_cone_iso_rot         = RR_X_90_CCW * gX_PT_els_in_cone_iso;            % rotate those points
gTH_PT_els_in_cone_iso_rot        = cartesian2spherical(gX_PT_els_in_cone_iso_rot); % convert back to spherical coordinates
lc_els_in_cone_iso                = local_coords_curved_3d_sph(GCOORD_SPH_rot_90_X,EL2NOD,els_in_cone_iso,gTH_PT_els_in_cone_iso_rot);

%===============================================================================================================================================================================================
% LOCAL COORDINATES FOR POINTS (gTH_PT)
%===============================================================================================================================================================================================
lc                                                 = zeros(3,size(gTH_PT,2));
lc(:,ismember(els,MESH.els_out_cone_no_cross_2pi)) = lc_els_out_cone_no_cross_2pi;
lc(:,ismember(els,MESH.els_out_cone_cross_2pi))    = lc_els_out_cone_cross_2pi;
lc(:,ismember(els,MESH.els_in_cone_no_iso))        = lc_els_in_cone_no_iso;
lc(:,ismember(els,MESH.els_in_cone_iso))           = lc_els_in_cone_iso;

end % END OF FUNCTION local_coords_3d_sph

% ##############################################################################################################################################################################################
%                               SUBFUNCTIONS
% ##############################################################################################################################################################################################

function lc = compute_lc(GCOORD_SPH,EL2NOD,els,gTH_PT)
% Usage: lc = compute_lc(GCOORD_SPH,EL2NOD,els,gTH_PT)
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

ns = length(els); % Number of points
th  = reshape(GCOORD_SPH(1,EL2NOD(1:4,els)),4,ns); % Coords of Vertex nodes
ph  = reshape(GCOORD_SPH(2,EL2NOD(1:4,els)),4,ns);
r   = reshape(GCOORD_SPH(3,EL2NOD(1:4,els)),4,ns);

thp = gTH_PT(1,:); % Global coords of points in each element
php = gTH_PT(2,:);
rp  = gTH_PT(3,:);

lc = zeros(3,length(els));

lc(1,:) =   (  (th(2,:)-thp)    .*(ph(3,:)-php)    .*(r(4,:) -rp)      + (th(3,:)-thp)    .*(ph(4,:)-php)    .*(r(2,:) -rp)      + (th(4,:)-thp)    .*(ph(2,:)-php)    .*(r(3,:) -rp)        ...
             - (r(2,:) -rp)     .*(ph(3,:)-php)    .*(th(4,:)-thp)     - (r(3,:) -rp)     .*(ph(4,:)-php)    .*(th(2,:)-thp)     - (r(4,:) -rp)     .*(ph(2,:)-php)    .*(th(3,:)-thp) )     ...
         ./ (  (th(2,:)-th(1,:)).*(ph(3,:)-ph(1,:)).*(r(4,:) -r(1,:))  + (th(3,:)-th(1,:)).*(ph(4,:)-ph(1,:)).*(r(2,:) -r(1,:))  + (th(4,:)-th(1,:)).*(ph(2,:)-ph(1,:)).*(r(3,:) -r(1,:))    ...
             - (r(2,:) -r(1,:)) .*(ph(3,:)-ph(1,:)).*(th(4,:)-th(1,:)) - (r(3,:) -r(1,:)) .*(ph(4,:)-ph(1,:)).*(th(2,:)-th(1,:)) - (r(4,:) -r(1,:)) .*(ph(2,:)-ph(1,:)).*(th(3,:)-th(1,:)) );

lc(2,:) =   (  (th(1,:)-thp)    .*(ph(3,:)-php)    .*(r(4,:) -rp)      + (th(3,:)-thp)    .*(ph(4,:)-php)    .*(r(1,:) -rp)      + (th(4,:)-thp)    .*(ph(1,:)-php)    .*(r(3,:) -rp)        ...
             - (r(1,:)-rp)      .*(ph(3,:)-php)    .*(th(4,:)-thp)     - (r(3,:) -rp)     .*(ph(4,:)-php)    .*(th(1,:)-thp)     - (r(4,:) -rp)     .*(ph(1,:)-php)    .*(th(3,:)-thp) )     ...
         ./ (  (th(1,:)-th(2,:)).*(ph(3,:)-ph(2,:)).*(r(4,:) -r(2,:))  + (th(3,:)-th(2,:)).*(ph(4,:)-ph(2,:)).*(r(1,:) -r(2,:))  + (th(4,:)-th(2,:)).*(ph(1,:)-ph(2,:)).*(r(3,:) -r(2,:))    ...
             - (r(1,:)-r(2,:))  .*(ph(3,:)-ph(2,:)).*(th(4,:)-th(2,:)) - (r(3,:) -r(2,:)) .*(ph(4,:)-ph(2,:)).*(th(1,:)-th(2,:)) - (r(4,:) -r(2,:)) .*(ph(1,:)-ph(2,:)).*(th(3,:)-th(2,:)) );
             
lc(3,:) =   (  (th(1,:)-thp)    .*(ph(2,:)-php)    .*(r(4,:) -rp)      + (th(2,:)-thp)    .*(ph(4,:)-php)    .*(r(1,:) -rp)      + (th(4,:)-thp)    .*(ph(1,:)-php)    .*(r(2,:) -rp)        ...
             - (r(1,:)-rp)      .*(ph(2,:)-php)    .*(th(4,:)-thp)     - (r(2,:) -rp)     .*(ph(4,:)-php)    .*(th(1,:)-thp)     - (r(4,:) -rp)     .*(ph(1,:)-php)    .*(th(2,:)-thp) )     ...
         ./ (  (th(1,:)-th(3,:)).*(ph(2,:)-ph(3,:)).*(r(4,:) -r(3,:))  + (th(2,:)-th(3,:)).*(ph(4,:)-ph(3,:)).*(r(1,:) -r(3,:))  + (th(4,:)-th(3,:)).*(ph(1,:)-ph(3,:)).*(r(2,:) -r(3,:))    ...
             - (r(1,:)-r(3,:))  .*(ph(2,:)-ph(3,:)).*(th(4,:)-th(3,:)) - (r(2,:) -r(3,:)) .*(ph(4,:)-ph(3,:)).*(th(1,:)-th(3,:)) - (r(4,:) -r(3,:)) .*(ph(1,:)-ph(3,:)).*(th(2,:)-th(3,:)) );

end % END OF SUBFUNCTION compute_lc