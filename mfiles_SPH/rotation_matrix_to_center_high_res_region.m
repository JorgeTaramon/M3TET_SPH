function RR_center = rotation_matrix_to_center_high_res_region(EP_lat,EP_lon,EP_angle)
% Usage: RR_center = rotation_matrix_to_center_high_res_region(EP_lat,EP_lon,EP_angle)
% 
% Purpose: 
%   Create a matrix to make a finite rotation from the rotation pole
%   (lat,lon) and finite rotation angle.
%
% Input:
%   EP_lat    : [scalar] : Euler Pole latitude (degrees)
%   EP_lon    : [scalar] : Euler Pole longitude (degrees)
%   EP_angle  : [scalar] : finite rotation angle (degrees)
%
% Output:
%   RR_center : [matrix] : Rotation matrix for Euler Pole
%
% JPM May 2011
% JMT May 2014: Fixed a bug related with the elements r_12 and r_21 from
%               the previous version. 
% JMT Feb 2017: Adapted to springmesh_3d

deg2rad           =  pi/180;
RR_center         = zeros(3,3);

s_EP_ang          = sin(deg2rad*EP_angle); % find sine of finite rotation angle
c_EP_ang          = cos(deg2rad*EP_angle); % find cos of finite rotation angle

vc_EP_ang         = 1 - c_EP_ang;

s_EP_lat          = sin(deg2rad*EP_lat); % find sine of Euler Pole latitude
c_EP_lat          = cos(deg2rad*EP_lat); % find cos of Euler Pole latitude
s_EP_lon          = sin(deg2rad*EP_lon); % find sine of Euler Pole longitude
c_EP_lon          = cos(deg2rad*EP_lon); % find cos of Euler Pole longitude

c_EP_lat_s_EP_lon = c_EP_lat .* sin(deg2rad*EP_lon); % c_EP_lat * s_EP_lon 
c_EP_lat_c_EP_lon = c_EP_lat .* cos(deg2rad*EP_lon); % c_EP_lat * c_EP_lon

RR_center(1,1,:)  = c_EP_ang + vc_EP_ang .* c_EP_lat_c_EP_lon.^2;
RR_center(1,2,:)  = vc_EP_ang .* c_EP_lat.^2 .* c_EP_lon .* s_EP_lon - s_EP_ang .* s_EP_lat;
RR_center(1,3,:)  = vc_EP_ang .* s_EP_lat .* c_EP_lat_c_EP_lon + s_EP_ang .* c_EP_lat_s_EP_lon;

RR_center(2,1,:)  = vc_EP_ang .* c_EP_lat.^2 .* c_EP_lon .* s_EP_lon + s_EP_ang .* s_EP_lat;
RR_center(2,2,:)  = c_EP_ang + vc_EP_ang .* c_EP_lat_s_EP_lon.^2;
RR_center(2,3,:)  = vc_EP_ang .* s_EP_lat .* c_EP_lat_s_EP_lon - s_EP_ang .* c_EP_lat_c_EP_lon;

RR_center(3,1,:)  = vc_EP_ang .* s_EP_lat .* c_EP_lat_c_EP_lon - s_EP_ang .* c_EP_lat_s_EP_lon;
RR_center(3,2,:)  = vc_EP_ang .* s_EP_lat .* c_EP_lat_s_EP_lon + s_EP_ang .* c_EP_lat_c_EP_lon;
RR_center(3,3,:)  = 1. - vc_EP_ang .* c_EP_lat.^2;

end % END OF FUNCTION rotation_matrix_to_center_high_res_region