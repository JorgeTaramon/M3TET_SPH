function [WS_els_in_cone_rot_X_90,WS_els_out_cone,WS_els_cross_2pi_rot_Z_180] = ...
    calc_tsearch2_ws_sph(GCOORD_SPH,EL2NOD,SETTINGS)
% Usage: [WS_els_in_cone_rot_X_90,WS_els_out_cone,WS_els_cross_2pi_rot_Z_180] = ...
%   calc_tsearch2_ws_sph(GCOORD_SPH,EL2NOD,SETTINGS)
% 
% Purpose:
%   Compute workspace for tsearch
%
% Input:
%   GCOORD_SPH : [matrix]    : spherical coordinates (theta,phi,r) 
%                             (in rad, rad and km) of the mesh
%   EL2NOD     : [matrix]    : connectivity matrix for GCOORD_SPH
%   SETTINGS   : [structure] : model parameters
%
% Output:
%   WS_els_in_cone_rot_X_90    : [structure] : workspace for points located
%                                              in elements within the cone
%   WS_els_out_cone            : [structure] : workspace for points located
%                                              in elements outside the cone
%   WS_els_cross_2pi_rot_Z_180 : [structure] : workspace for points located
%                                              in elements crossing phi =
%                                              2pi
%
% JMT Oct 2016

%=========================================================================================================================================================================================================
% COMPUTE ELEMENTS IN RELATION WITH THE CONE AND CROSSING PHI = 2PI
%=========================================================================================================================================================================================================
if iscell(EL2NOD)
    EL2NOD = EL2NOD{1};
end
[els_in_cone,~,~,els_out_cone,els_in_cone_iso] = check_els_in_cone(GCOORD_SPH,EL2NOD,SETTINGS.theta_cone);
[els_cross_2pi,~]         = check_phi(GCOORD_SPH,EL2NOD(:,els_out_cone));
els_out_cone_cross_2pi    = els_out_cone(els_cross_2pi);
els_out_cone_no_cross_2pi = els_out_cone;
els_out_cone_no_cross_2pi(ismember(els_out_cone_no_cross_2pi,els_out_cone_cross_2pi)) = [];
els_in_cone_no_iso        = els_in_cone;
els_in_cone_no_iso(ismember(els_in_cone_no_iso,els_in_cone_iso)) = [];
GCOORD                    = spherical2cartesian(GCOORD_SPH);

%=========================================================================================================================================================================================================
% COMPUTE INITIAL GUESS FOR POINTS TO BE LOCATED
%=========================================================================================================================================================================================================
% 1.Compute initial guess for points inside elements that are outise the cone and not crossing phi = 2pi
th_PT_els_out_cone_no_cross_2pi   = mean( reshape(GCOORD_SPH(1,EL2NOD(1:4,els_out_cone_no_cross_2pi)),4,size(els_out_cone_no_cross_2pi,2)) ,1);
ph_PT_els_out_cone_no_cross_2pi   = mean( reshape(GCOORD_SPH(2,EL2NOD(1:4,els_out_cone_no_cross_2pi)),4,size(els_out_cone_no_cross_2pi,2)) ,1);
r_PT_els_out_cone_no_cross_2pi    = mean( reshape(GCOORD_SPH(3,EL2NOD(1:4,els_out_cone_no_cross_2pi)),4,size(els_out_cone_no_cross_2pi,2)) ,1);
gTH_PT_els_out_cone_no_cross_2pi  = [th_PT_els_out_cone_no_cross_2pi; ph_PT_els_out_cone_no_cross_2pi;r_PT_els_out_cone_no_cross_2pi];
clear th_PT_els_out_cone_no_cross_2pi ph_PT_els_out_cone_no_cross_2pi r_PT_els_out_cone_no_cross_2pi

% 2.Compute initial guess for points inside elements that are outise the cone and crossing phi = 2pi
RR_Z_180_CCW                      = [-1  0  0 ; ...
                                      0 -1  0 ; ...
                                      0  0  1 ];
GCOORD_rot_Z_180                  = RR_Z_180_CCW * GCOORD;
GCOORD_SPH_rot_Z_180              = cartesian2spherical(GCOORD_rot_Z_180);
th_PT_els_out_cone_cross_2pi_rot  = mean( reshape(GCOORD_SPH_rot_Z_180(1,EL2NOD(1:4,els_out_cone_cross_2pi)),4,size(els_out_cone_cross_2pi,2)) ,1);
ph_PT_els_out_cone_cross_2pi_rot  = mean( reshape(GCOORD_SPH_rot_Z_180(2,EL2NOD(1:4,els_out_cone_cross_2pi)),4,size(els_out_cone_cross_2pi,2)) ,1);
r_PT_els_out_cone_cross_2pi_rot   = mean( reshape(GCOORD_SPH_rot_Z_180(3,EL2NOD(1:4,els_out_cone_cross_2pi)),4,size(els_out_cone_cross_2pi,2)) ,1);
gTH_PT_els_out_cone_cross_2pi_rot = [th_PT_els_out_cone_cross_2pi_rot; ph_PT_els_out_cone_cross_2pi_rot;r_PT_els_out_cone_cross_2pi_rot];
gX_PT_els_out_cone_cross_2pi_rot  = spherical2cartesian(gTH_PT_els_out_cone_cross_2pi_rot);
RR_Z_180_CW                       = [-1  0  0 ; ...
                                      0 -1  0 ; ...
                                      0  0  1 ];
gX_PT_els_out_cone_cross_2pi      = RR_Z_180_CW * gX_PT_els_out_cone_cross_2pi_rot;
gTH_PT_els_out_cone_cross_2pi     = cartesian2spherical(gX_PT_els_out_cone_cross_2pi);
clear th_PT_els_out_cone_cross_2pi_rot ph_PT_els_out_cone_cross_2pi_rot r_PT_els_out_cone_cross_2pi_rot gX_PT_els_out_cone_cross_2pi_rot

% 3.Compute initial guess for points inside elements that are within the cone and are not isoparametric
RR_X_90_CCW                       = [ 1  0  0 ; ...
                                      0  0 -1 ; ...
                                      0  1  0 ];
GCOORD_rot_X_90                   = RR_X_90_CCW * GCOORD;
GCOORD_SPH_rot_X_90               = cartesian2spherical(GCOORD_rot_X_90);
th_PT_els_in_cone_no_iso_rot      = mean( reshape(GCOORD_SPH_rot_X_90(1,EL2NOD(1:4,els_in_cone_no_iso)),4,size(els_in_cone_no_iso,2)) ,1);
ph_PT_els_in_cone_no_iso_rot      = mean( reshape(GCOORD_SPH_rot_X_90(2,EL2NOD(1:4,els_in_cone_no_iso)),4,size(els_in_cone_no_iso,2)) ,1);
r_PT_els_in_cone_no_iso_rot       = mean( reshape(GCOORD_SPH_rot_X_90(3,EL2NOD(1:4,els_in_cone_no_iso)),4,size(els_in_cone_no_iso,2)) ,1);
gTH_PT_els_in_cone_no_iso_rot     = [th_PT_els_in_cone_no_iso_rot; ph_PT_els_in_cone_no_iso_rot;r_PT_els_in_cone_no_iso_rot];
gX_PT_els_in_cone_no_iso_rot      = spherical2cartesian(gTH_PT_els_in_cone_no_iso_rot);
RR_X_90_CW                        = [ 1  0  0 ; ...
                                      0  0  1 ; ...
                                      0 -1  0 ];
gX_PT_els_in_cone_no_iso          = RR_X_90_CW * gX_PT_els_in_cone_no_iso_rot;
gTH_PT_els_in_cone_no_iso         = cartesian2spherical(gX_PT_els_in_cone_no_iso);
clear th_PT_els_in_cone_no_iso_rot ph_PT_els_in_cone_no_iso_rot r_PT_els_in_cone_no_iso_rot gX_PT_els_in_cone_no_iso_rot

% 4.Compute initial guess for points inside elements that are within the cone and are isoparametric
RR_X_90_CCW                       = [ 1  0  0 ; ...
                                      0  0 -1 ; ...
                                      0  1  0 ];
GCOORD_rot_X_90                   = RR_X_90_CCW * GCOORD;
GCOORD_SPH_rot_X_90               = cartesian2spherical(GCOORD_rot_X_90);
th_PT_els_in_cone_iso_rot         = mean( reshape(GCOORD_SPH_rot_X_90(1,EL2NOD(1:4,els_in_cone_iso)),4,size(els_in_cone_iso,2)) ,1);
ph_PT_els_in_cone_iso_rot         = mean( reshape(GCOORD_SPH_rot_X_90(2,EL2NOD(1:4,els_in_cone_iso)),4,size(els_in_cone_iso,2)) ,1);
r_PT_els_in_cone_iso_rot          = mean( reshape(GCOORD_SPH_rot_X_90(3,EL2NOD(1:4,els_in_cone_iso)),4,size(els_in_cone_iso,2)) ,1);
gTH_PT_els_in_cone_iso_rot        = [th_PT_els_in_cone_iso_rot; ph_PT_els_in_cone_iso_rot;r_PT_els_in_cone_iso_rot];
gX_PT_els_in_cone_iso_rot         = spherical2cartesian(gTH_PT_els_in_cone_iso_rot);
RR_X_90_CW                        = [ 1  0  0 ; ...
                                      0  0  1 ; ...
                                      0 -1  0 ];
gX_PT_els_in_cone_iso             = RR_X_90_CW * gX_PT_els_in_cone_iso_rot;
gTH_PT_els_in_cone_iso            = cartesian2spherical(gX_PT_els_in_cone_iso);
clear th_PT_els_in_cone_iso_rot ph_PT_els_in_cone_iso_rot r_PT_els_in_cone_iso_rot gX_PT_els_in_cone_iso_rot

% 5.Compute initial guess for points to be located
gTH_PT = [gTH_PT_els_out_cone_no_cross_2pi ...
          gTH_PT_els_out_cone_cross_2pi    ...
          gTH_PT_els_in_cone_no_iso        ...
          gTH_PT_els_in_cone_iso];         % points to be located
gTH_PT = [gTH_PT ...
          gTH_PT_els_out_cone_cross_2pi(:, [1 4 7])+0.01 ....
          gTH_PT_els_out_cone_no_cross_2pi(:,[4 9 2])+0.01 ...
          gTH_PT_els_out_cone_no_cross_2pi(:,[2 6 1])+0.01 ...
          gTH_PT_els_in_cone_no_iso(:,[9 7 4])+0.01];
gX_PT  = spherical2cartesian(gTH_PT);
% clear gTH_PT_els_out_cone_no_cross_2pi gTH_PT_els_out_cone_cross_2pi gTH_PT_els_in_cone_no_iso gTH_PT_els_in_cone_iso

if size(gTH_PT,2) < 100
    ligthWhi  = [1 1 1];       % colour for elements outside the cone and not crossing phi = 2pi
    gX_PT_els_out_cone_no_cross_2pi = spherical2cartesian(gTH_PT_els_out_cone_no_cross_2pi);
    plot_points_and_elements(GCOORD,GCOORD_SPH,EL2NOD,gX_PT_els_out_cone_no_cross_2pi,gTH_PT_els_out_cone_no_cross_2pi,els_out_cone_no_cross_2pi,SETTINGS.theta_cone,ligthWhi,95)
    subplot(1,2,1)
    title('Points to be located inside elements outside the cone and not crossing \phi = 2\pi in Cartersian coordinates')
    subplot(1,2,2)
    title('Points to be located inside elements outside the cone and not crossing \phi = 2\pi in spherical coordinates')
    
    lightYell = 0.8*[1 1 0];   % colour for elements outside the cone and crossing phi = 2pi
    plot_points_and_elements(GCOORD,GCOORD_SPH_rot_Z_180,EL2NOD,gX_PT_els_out_cone_cross_2pi,gTH_PT_els_out_cone_cross_2pi_rot,els_out_cone_cross_2pi,SETTINGS.theta_cone,lightYell,96)
    subplot(1,2,1)
    title('Points to be located inside elements outside the cone and crossing \phi = 2\pi in Cartersian coordinates')
    subplot(1,2,2)
    title('Points to be located inside elements outside the cone and crossing \phi = 2\pi in 180° rotated spherical frame')
    
    lightMag  = 0.8*[1 0.7 1]; % colour for the elements within the cone and not isoparametric
    plot_points_and_elements(GCOORD,GCOORD_SPH_rot_X_90,EL2NOD,gX_PT_els_in_cone_no_iso,gTH_PT_els_in_cone_no_iso_rot,els_in_cone_no_iso,SETTINGS.theta_cone,lightMag,97)
    subplot(1,2,1)
    title('Points to be located inside elements within the cone and not isoparametric in Cartersian coordinates')
    subplot(1,2,2)
    title('Points to be located inside elements within the cone and not isoparametric in 90° rotated spherical frame')
    
    lightCyan = 0.8*[0 1 1];   % colour for the elements within the cone and isoparametric
    plot_points_and_elements(GCOORD,GCOORD_SPH_rot_X_90,EL2NOD,gX_PT_els_in_cone_iso,gTH_PT_els_in_cone_iso_rot,els_in_cone_iso,SETTINGS.theta_cone,lightCyan,98)
    subplot(1,2,1)
    title('Points to be located inside elements within the cone and isoparametric in Cartersian coordinates')
    subplot(1,2,2)
    title('Points to be located inside elements within the cone and isoparametric in 90° rotated spherical frame')
end

%=========================================================================================================================================================================================================
% LOCATE POINTS INSIDE ELEMENTS 
%=========================================================================================================================================================================================================
calc_element_numbers = 'yes';
if strcmp(calc_element_numbers,'no')
    % only search for a small number of points because only the
    % tsearch2-workspaces are of interest here
    nPT    = 100;
    gTH_PT = gTH_PT(:,1:nPT);
    gX_PT  = gX_PT(:,1:nPT);
end

theta_cone     = SETTINGS.theta_cone*pi/180; % theta cone in rad
theta_north    = theta_cone;                 % theta angle from the +Z axis to the generatrix
theta_south    = pi-theta_cone;              % theta angle from the -Z axis to the generatrix
nodes_in_cone  = gTH_PT(1,:) <= theta_north | gTH_PT(1,:) >= theta_south; % boolean vector for nodes in the cone
nnod           = max(max(EL2NOD(1:4,:)));
els            = zeros(size(gTH_PT,2),3);

% LOCATE POINTS INSIDE ELEMENTS WITHIN THE CONE
gX_PT_rot_X_90          = RR_X_90_CCW * gX_PT;
gTH_PT_rot_X_90         = cartesian2spherical(gX_PT_rot_X_90);
gTH_PT_in_cone_rot_X_90 = gTH_PT_rot_X_90(:,nodes_in_cone);
[els_in_cone_rot_X_90,WS_els_in_cone_rot_X_90,~] = ...
    tsearch2(GCOORD_SPH_rot_X_90(:,1:nnod),EL2NOD([1 2 4 3],:),gTH_PT_in_cone_rot_X_90,[],[]);
els(nodes_in_cone,1)  = els_in_cone_rot_X_90';

% LOCATE POINTS INSIDE ELEMENTS OUTSIDE THE CONE
gTH_PT_out_cone = gTH_PT(:,~nodes_in_cone);
[els_out_cone,WS_els_out_cone,~] = ...
    tsearch2(GCOORD_SPH(:,1:nnod),EL2NOD([1 2 4 3],:),gTH_PT_out_cone,[],[]);
els(~nodes_in_cone,2) = els_out_cone';

% LOCATE POINTS INSIDE ELEMENTS CROSSING PHI = 2PI
nodes_cross_2pi  = sum(els,2) == 0;
gX_PT_rot_Z_180  = RR_Z_180_CCW * gX_PT;
gTH_PT_rot_Z_180 = cartesian2spherical(gX_PT_rot_Z_180);
gTH_PT_cross_2pi_rot_Z_180 = gTH_PT_rot_Z_180(:,nodes_cross_2pi);
[els_cross_2pi_rot_Z_180,WS_els_cross_2pi_rot_Z_180,~] = ...
    tsearch2(GCOORD_SPH_rot_Z_180(:,1:nnod),EL2NOD([1 2 4 3],:),gTH_PT_cross_2pi_rot_Z_180,[],[]);
els(nodes_cross_2pi,3) = els_cross_2pi_rot_Z_180';
els = uint32(sum(els,2));

end

% ########################################################################################################################################################################################################
%                                           SUB-FUNCTIONS
% ########################################################################################################################################################################################################

function plot_points_and_elements(GCOORD,GCOORD_SPH,EL2NOD,gX_PT,gTH_PT,els,theta_cone,color,Fig_num)

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
tetramesh(EL2NOD(1:4,els)',GCOORD','FaceColor',color,'FaceAlpha',0.3)
scatter3(gX_PT(1,:),gX_PT(2,:),gX_PT(3,:),'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
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
tetramesh(EL2NOD(1:4,els)',GCOORD_SPH_plot','FaceColor',color,'FaceAlpha',0.3)
scatter3(gTH_PT(1,:),gTH_PT(2,:),gTH_PT(3,:)/1e3,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
xlabel('\theta (rad)')
ylabel('\phi (rad)')
zlabel('r (10^3 km)')

end % END OF SUBFUNCTION plot_points_and_elements