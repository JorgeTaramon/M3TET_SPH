function [VAR,NUMSCALE] = init_fields_sph(MESH,SETTINGS,PHYSICS,NUMSCALE)
% Usage: [VAR,NUMSCALE] = init_fields_sph(MESH,SETTINGS,PHYSICS,NUMSCALE)
% 
% Purpose: Initialize major variable fields and store in structure "VAR".
%
% Input:
%   MESH     : [structure] : FE mesh parameters
%   SETTINGS : [structure] : model parameters
%   PHYSICS  : [structure] : physical properties
%   NUMSCALE : [structure] : numerical scaling parameters
%
% Output:
%   VAR      : [structure] : major variable fields, each is a vector
%   NUMSCALE : [structure] : numerical scaling parameters
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Jan 2011
% JMT Feb 2017: Added initial temperature field for cratons in South
%               Atlantic region

n   = MESH.nnod;   % number of temperature unknowns
nP  = MESH.nVnod;  % number of pressure nodes

% Allocate memory for all variables
VAR.Ux      = zeros(n,1);  % x-velocity component
VAR.Uy      = zeros(n,1);  % y-velocity component
VAR.Uz      = zeros(n,1);  % z-velocity component
if strcmp(SETTINGS.plate_model,'WJM')
    VAR.Ur  = zeros(n,1); % r-velocity component
    VAR.Ue  = zeros(n,1); % east-velocity component
    VAR.Un  = zeros(n,1); % north-velocity component
elseif sum(strcmp(SETTINGS.plate_model,{'GPlates','Debug'}))
    VAR.Uth = zeros(n,1); % theta-velocity component
    VAR.Uph = zeros(n,1); % phi-velocity component
    VAR.Ur  = zeros(n,1); % r-velocity component
end
VAR.P       = zeros(nP,1); % pressure unknowns (1 per vertex node)
VAR.T       = zeros(n,1);  % temperature unknowns (1 per node)
VAR.Dens    = PHYSICS.Dens*ones(n,1); % density (1 per node, depending on T)
VAR.Visc    = PHYSICS.Visc*ones(n,1); % viscosity (1 per node, depending on T)
VAR.ref     = MESH.points_ref;

% Load coordinates
x                    = MESH.GCOORD(1,:)';
y                    = MESH.GCOORD(2,:)';
z                    = MESH.GCOORD(3,:)';
r                    = MESH.GCOORD_SPH(3,:);
r(MESH.PointID==301) = MESH.r_cmb; 
r(MESH.PointID==306) = MESH.r_surf;

% Load top, bottom and mantle temperatures
Ttop = PHYSICS.Ttop;
Tbot = PHYSICS.Tbot;
Tm   = PHYSICS.Tm;

%=================================================================================================================================
% DEFINE INITIAL TEMPERATURE FIELD (Comment/uncomment the following options to choose a initial temperature field) 
%=================================================================================================================================

%     %======================================================================
%     % LINEAR TEMPERATURE GRADIENT
%     m     = (Ttop - Tbot)/(MESH.r_surf-MESH.r_cmb);
%     b     = Ttop - m*MESH.r_surf;
%     VAR.T = m*r + b;
%     %======================================================================
    
%     %======================================================================
%     % UNIFORM TEMPERATURE EVERYWHERE
%     VAR.T(:) = Tm;
%     %======================================================================
    
%     %======================================================================
%     % ADD HOT LAYER AT BOTTOM
%      VAR.T(r<MESH.r_cmb+300) = Tbot;
%     %======================================================================
    
%     %======================================================================
%     % ADD COLD LAYER AT TOP
%     VAR.T(r>MESH.r_surf-200) = Ttop;
%     %======================================================================
    
%     %======================================================================
%     % ADD 2 HOT BLOBS
%     % Blob 1
%     xc          = 0;
%     yc          = 4921;
%     zc          = 0;
%     r           = 500;
%     ind1        = sqrt( (x-xc).^2 + (y-yc).^2 + (z-zc).^2 ) <= r;
%     VAR.T(ind1) = VAR.T(ind1) + 150;
%     % Blob 2
%     xc          = 0;
%     yc          = -4921;
%     zc          = 0;
%     r           = 500;
%     ind2        = sqrt( (x-xc).^2 + (y-yc).^2 + (z-zc).^2 ) <= r;
%     VAR.T(ind2) = VAR.T(ind2) + 150;
%     %======================================================================
    
%     %======================================================================
%     % ADD HOT BLOB AT CMB
%     xc          = 0;
%     yc          = 3471;
%     zc          = 0;
%     r           = 800;
%     ind1        = sqrt( (x-xc).^2 + (y-yc).^2 + (z-zc).^2 ) <= r;
%     VAR.T(ind1) = VAR.T(ind1) + 300;
%     %======================================================================
    
%     %===============================================================================================================
%     % TEMPERATURE PROFILE FOR CONTINENTAL LITHOSPHERE EVERYWHERE (age ~200Myr)(USING HSCM)
%     ylim                       = 1e3 * [3471 6371] ./ NUMSCALE.L0;
%     Myr                        = 1e6*365.25*24*3600;
%     cont_age_SI                = 200 * Myr;  % target continental isotherm
%     dt_SI                      = 10 * Myr; % time step
%     nstep_continental_isotherm = cont_age_SI/dt_SI; % number of time steps
%     time_SI                    = 0; % initial time
%     [~,p]                      = sort(r*NUMSCALE.L0);
%     figure(57); clf
%     for i = 1:nstep_continental_isotherm
%         T       = solution_halfspace_cool(r*NUMSCALE.L0,MESH.r_surf*NUMSCALE.L0,PHYSICS.kappa_SI,time_SI,Ttop,Tbot);
%         time_SI = time_SI + dt_SI;
%     end
%     plot(T(p),r(p),'r.-')
%     title('Isotherm for continental crust')
%     set(gca,'YLim',ylim);
%     xlabel('T (C)');
%     ylabel(sprintf('Distance (%s)',NUMSCALE.unit_L));
%     grid on
%     VAR.T = T';
%     %===============================================================================================================

    %=============================================================================================================================
    % INITIAL TEMPERATURE FOR SOUTH ATLANTIC OPENING (130Myr-100Myr) USING HSCM
    Myr                    = NUMSCALE.t0;
    if strcmp(SETTINGS.cratons_contours,'3')
        ocean_age_SI       =  80 * Myr; % ocean age
        cont_age_SI        = 100 * Myr; % continental age
        craton_age_SI_bnd  = [200 275 350] * Myr; % contour cratons ages for smoothing the shape
    elseif strcmp(SETTINGS.cratons_contours,'8') || strcmp(SETTINGS.cratons_contours,'8M')
        ocean_age_SI       = 100 * Myr; % ocean age
        cont_age_SI        = 100 * Myr; % continental age
        craton_age_SI_bnd  = [150 200 250 275 300 325 340 350] * Myr; % contour cratons ages for smoothing the shape
    end

    % - Outside the refined region -> Temperature profile for oceaninc lithosphere
    % -----------------------------------------------------------------------------------------
        nodes_ocean = 1:size(MESH.GCOORD,2);
        figure(56); clf
        subplot(1,3,1)
        T_ocean = HSCM(ocean_age_SI,nodes_ocean,NUMSCALE,MESH,PHYSICS);
        title('Isotherm for oceanic lithosphere')
        VAR.T(nodes_ocean) = T_ocean;
        figure(57); clf
        subplot(1,3,1)
        plot_density(T_ocean',nodes_ocean,MESH,PHYSICS,NUMSCALE);
        title('Density for oceanic lithosphere')
        figure(58); clf
        subplot(1,3,1)
        plot_viscosity(T_ocean',nodes_ocean,MESH,PHYSICS,NUMSCALE);
        title('Viscosity for oceanic lithosphere')
    
    % - Inside the refined region -> Temperature profile for continental lithosphere
    % --------------------------------------------------------------------------------------------
        deg2rad          = pi/180;
        theta0           = MESH.theta0;                 % colatitude (degrees) of the point around which refined zones are defined
        phi0             = MESH.phi0;                   % longitude (degrees) of the point around which refined zones are defined
        w_ref            = MESH.w_ref;                  % width of the refined zone (km)
        w_ref_deg        = w_ref/(deg2rad*MESH.r_surf); % width of refined zone in degrees (North-South)
        theta_ref_n      = theta0-w_ref_deg/2;          % colatitude of the northern boundary in the refined zone
        theta_ref_s      = theta0+w_ref_deg/2;          % colatitude of the southern boundary in the refined zone
        l_ref            = MESH.l_ref;                  % length of the refined zone (km)
        l_ref_deg        = l_ref/(deg2rad*MESH.r_surf); % length of refined zone in degrees (East-West)
        phi_ref_e        = phi0+l_ref_deg/2;            % longitude of the eastern boundary in the refined zone
        phi_ref_w        = phi0-l_ref_deg/2;            % longitude of the westren boundary in the refined zone
        if SETTINGS.del_lith_thickness == 0
            lat_lon_refined  = [90-theta_ref_n 90-theta_ref_n 90-theta_ref_s 90-theta_ref_s; ...
                                phi_ref_w      phi_ref_e      phi_ref_e      phi_ref_w];
            nodes_cont       = nodes_inside_refined_region(MESH,lat_lon_refined);
            figure(56)
            subplot(1,3,2)
            T_cont = HSCM(cont_age_SI,nodes_cont,NUMSCALE,MESH,PHYSICS);
            title('Isotherm for continental lithosphere')
            VAR.T(nodes_cont) = T_cont;
            figure(57)
            subplot(1,3,2)
            plot_density(T_cont',nodes_ocean(nodes_cont),MESH,PHYSICS,NUMSCALE);
            title('Density for continental lithosphere')
            figure(58)
            subplot(1,3,2)
            plot_viscosity(T_cont',nodes_ocean(nodes_cont),MESH,PHYSICS,NUMSCALE);
            title('Viscosity for continental lithosphere')
        else
            % Compute colatitude of the plume in the rotated frame (where high-res is centered in the equator) 
            colat_lon_plume_center                 = PHYSICS.PLUME.colat_lon_plume_center; % center of the plume (colatitude,longitude) [degrees]
            theta_plume_center                     = colat_lon_plume_center(1)*pi/180; % colatitude in radians
            phi_plume_center                       = colat_lon_plume_center(2)*pi/180; % longitude in radians
            phi_plume_center(phi_plume_center < 0) = phi_plume_center(phi_plume_center < 0) + 2*pi; % values between 0 and 2pi
            GCOORD_plume_center                    = spherical2cartesian([theta_plume_center; phi_plume_center; MESH.r_surf]);
            GCOORD_plume_center_rot                = (MESH.RR2 * MESH.RR1)' * GCOORD_plume_center; % rotate coordinates of the plume center from GPlates frame
            GCOORD_SPH_plume_center_rot            = cartesian2spherical(GCOORD_plume_center_rot);
            theta_plume_rot                        = GCOORD_SPH_plume_center_rot(1,1)*180/pi;
            
            % Compute initial temperature profile from plume location towards the North  
            lat_lon_refined_north = [90-theta_ref_n 90-theta_ref_n 90-theta_plume_rot 90-theta_plume_rot; ...
                                     phi_ref_w      phi_ref_e      phi_ref_e          phi_ref_w         ];
            nodes_cont_north      = nodes_inside_refined_region(MESH,lat_lon_refined_north);
            figure(56)
            subplot(1,3,2)
            T_cont_north = HSCM(cont_age_SI,nodes_cont_north,NUMSCALE,MESH,PHYSICS);
            title('Isotherm for continental lithosphere')
            VAR.T(nodes_cont_north) = T_cont_north;
            figure(57)
            subplot(1,3,2)
            plot_density(T_cont_north',nodes_ocean(nodes_cont_north),MESH,PHYSICS,NUMSCALE);
            title('Density for continental lithosphere')
            figure(58)
            subplot(1,3,2)
            plot_viscosity(T_cont_north',nodes_ocean(nodes_cont_north),MESH,PHYSICS,NUMSCALE);
            title('Viscosity for continental lithosphere')
            
            % Compute initial temperature profile from plume location towards the South
            lat_lon_refined_south = [90-theta_plume_rot 90-theta_plume_rot 90-theta_ref_s 90-theta_ref_s; ...
                                     phi_ref_w          phi_ref_e          phi_ref_e      phi_ref_w     ];
            nodes_cont_south      = nodes_inside_refined_region(MESH,lat_lon_refined_south);
            
            % compute continental age in the shouthern end of the high-res region 
            del_lith_thickness = SETTINGS.del_lith_thickness;
            z_cont_north       = 2.32 * sqrt(PHYSICS.kappa_SI * cont_age_SI)/1e3; % thickness of continental lithosphere (in km) 
                                                                                  % in the northern part of the high-res region
            z_cont_south       = z_cont_north + del_lith_thickness; % thickness of continental lithosphere (in km) in the southern part of the high-res region
            cont_age_south     = z_cont_south^2 * 1e6/(2.32^2 * PHYSICS.kappa_SI * Myr);
            
            % compute slope for the change of age as a function of colatitude (Myr/rad) 
            m = (cont_age_south - cont_age_SI/Myr)/((theta_ref_s - theta_plume_rot)*pi/180);
            b = cont_age_SI/Myr - (m * theta_plume_rot*pi/180);
            
            % compute continental age for nodes 
            theta_nodes_south = MESH.GCOORD_SPH(1,nodes_cont_south);
            cont_age_south    = m*theta_nodes_south + b;
            cont_age_south_SI = cont_age_south * Myr;
            
            figure(56)
            subplot(1,3,2)
            T_cont_south = HSCM(cont_age_south_SI,nodes_cont_south,NUMSCALE,MESH,PHYSICS);
            title('Isotherm for continental lithosphere')
            VAR.T(nodes_cont_south) = T_cont_south;
            figure(57)
            subplot(1,3,2)
            plot_density(T_cont_south',nodes_ocean(nodes_cont_south),MESH,PHYSICS,NUMSCALE);
            title('Density for continental lithosphere')
            figure(58)
            subplot(1,3,2)
            plot_viscosity(T_cont_south',nodes_ocean(nodes_cont_south),MESH,PHYSICS,NUMSCALE);
            title('Viscosity for continental lithosphere')
        end
        
    
    % - Inside the cratons -> Temperature profile for cratons
    % ----------------------------------------------------------------------
        n_contours    = size(PHYSICS.CRATONS,2);
        if strcmp(SETTINGS.cratons_contours,'3') || strcmp(SETTINGS.cratons_contours,'8')
            n_cratons = 5;
        elseif strcmp(SETTINGS.cratons_contours,'8M')
            n_cratons = 3;
        end
        n_countours_per_craton           = n_contours/n_cratons;
        nodes_inside_craton(1:n_cratons) = struct('nodes',0);
        figure(56)
        subplot(1,3,3)
        for i = 1:n_cratons
            for j = 1:n_countours_per_craton
                nodes_inside_craton(n_countours_per_craton*(i-1)+j).nodes = inside_craton(MESH,PHYSICS.CRATONS(n_countours_per_craton*(i-1)+j).latlon');
                nodes_craton                                   = 1:size(MESH.GCOORD,2);
                VAR.T(nodes_inside_craton(n_countours_per_craton*(i-1)+j).nodes) = ...
                    HSCM(craton_age_SI_bnd(j),nodes_craton(nodes_inside_craton(n_countours_per_craton*(i-1)+j).nodes),NUMSCALE,MESH,PHYSICS);
            end
        end
        title('Isotherm for craton')
        figure(57)
        subplot(1,3,3)
        plot_density(VAR.T(nodes_inside_craton(n_countours_per_craton*(i-1)+j).nodes),nodes_ocean(nodes_inside_craton(n_countours_per_craton*(i-1)+j).nodes),MESH,PHYSICS,NUMSCALE);
        title('Density for craton')
        figure(58)
        subplot(1,3,3)
        plot_viscosity(VAR.T(nodes_inside_craton(n_countours_per_craton*(i-1)+j).nodes),nodes_ocean(nodes_inside_craton(n_countours_per_craton*(i-1)+j).nodes),MESH,PHYSICS,NUMSCALE);
        title('Viscosity for craton')
    %=============================================================================================================================

%====================================================================================================================================================
% DEFINE NUMSCALE.Dens0 
%====================================================================================================================================================

    %=================================================================================================================================
    % TEMPERATURE-DEPENDENT DENSITY AT THE BOTTOM OF REFINED REGION FOLLOWING THE TEMPERATURE PROFILE FOR A CONTINENTAL LITHOSPHERE
    cont_age_SI         = 100 * Myr;  % target continental isotherm
    T_bot_refied_region = ...
        solution_halfspace_cool((MESH.r_surf-MESH.d_ref)*NUMSCALE.L0,MESH.r_surf*NUMSCALE.L0,PHYSICS.kappa_SI,cont_age_SI,Ttop,Tbot);
    NUMSCALE.Dens0      = PHYSICS.Dens * (1 - PHYSICS.alpha*(T_bot_refied_region-PHYSICS.T_ref));
    %================================================================================================================================

    %====================================================================================================================================================
    % TEMPERATURE-DEPENDENT DENSITY FOR BARYCENTERS OF ELEMENTS INSIDE THE REFINED REGION FOLLOWING THE TEMPERATURE PROFILE FOR A CONTINENTAL LITHOSPHERE
    GCOORD      = MESH.GCOORD;
    EL2NOD      = MESH.EL2NOD{1};
    els_ref     = MESH.els_ref;
    x_bary      = (GCOORD(1,EL2NOD(1,els_ref)) + GCOORD(1,EL2NOD(2,els_ref)) + GCOORD(1,EL2NOD(3,els_ref)) + GCOORD(1,EL2NOD(4,els_ref)))/4;
    y_bary      = (GCOORD(2,EL2NOD(1,els_ref)) + GCOORD(2,EL2NOD(2,els_ref)) + GCOORD(2,EL2NOD(3,els_ref)) + GCOORD(2,EL2NOD(4,els_ref)))/4;
    z_bary      = (GCOORD(3,EL2NOD(1,els_ref)) + GCOORD(3,EL2NOD(2,els_ref)) + GCOORD(3,EL2NOD(3,els_ref)) + GCOORD(3,EL2NOD(4,els_ref)))/4;
    r_bary      = sqrt(x_bary.^2+y_bary.^2+z_bary.^2);
    cont_age_SI = 100 * Myr;  % target continental isotherm
	T_ref_bary  = solution_halfspace_cool(r_bary*NUMSCALE.L0,MESH.r_surf*NUMSCALE.L0,PHYSICS.kappa_SI,cont_age_SI,Ttop,Tbot);
    NUMSCALE.Dens0_T_bary = PHYSICS.Dens * (1 - PHYSICS.alpha*(T_ref_bary' - PHYSICS.T_ref));
    figure(59); clf
    plot(NUMSCALE.Dens0_T_bary,r_bary,'r.')
    title('Density for continental crust')
    set(gca,'YLim',ylim);
    xlabel('\rho (kg/m^3)');
    ylabel(sprintf('Distance (%s)',NUMSCALE.unit_L));
    grid on
    %====================================================================================================================================================

end % END OF FUNCTION init_fields_sph

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function T = HSCM(age_SI,nodes,NUMSCALE,MESH,PHYSICS)
r       = MESH.GCOORD_SPH(3,nodes);
[~,p]   = sort(r*NUMSCALE.L0);
ylim    = 1e3 * [MESH.r_cmb MESH.r_surf] ./ NUMSCALE.L0;
T       = solution_halfspace_cool(r*NUMSCALE.L0,MESH.r_surf*NUMSCALE.L0,PHYSICS.kappa_SI,age_SI,PHYSICS.Ttop,PHYSICS.Tbot);
plot(T(p),r(p),'r.-')
set(gca,'YLim',ylim);
xlabel('T (C)');
ylabel(sprintf('Distance (%s)',NUMSCALE.unit_L));
grid on

end

% #########################################################################

function T = solution_halfspace_cool(z,z_T_top,kappa,time,T_top,T_bot)

scale_time = 2*sqrt(kappa*time);
if scale_time<=0
    T             = T_bot * ones(size(z));
    T(z==z_T_top) = T_top;
else
    T = T_top + (T_bot-T_top)*erf(abs(z_T_top-z)./scale_time);
end

end % END OF SUBFUNCTION solution_halfspace_cool

% #########################################################################

function plot_viscosity(T,nodes,MESH,PHYSICS,NUMSCALE)

% Reference values
% =================
vALL  = ones(size(nodes,2),1);
R     = 8.314472; % J/mol K; univ. gas constant
Visc0 = PHYSICS.Visc;
Dens0 = PHYSICS.Dens;
T0    = PHYSICS.T_ref;
r     = MESH.GCOORD_SPH(3,nodes);
[~,p] = sort(r*NUMSCALE.L0);
depth = MESH.r_surf - r;

% Temperature effect on viscosity
% ===============================
T0_K  = T0 + 273;   % T0 in Kelvin
Ea    = PHYSICS.Ea; % J/mol; activation energy
if Ea>0
    if Ea<1000
        error(' Activation energy PHYSICS.Ea seems to be in kJ/mol. Convert to J/mol.');
    end
    T_K  = T + 273;   % T in Kelvin
    %     T_K  = T_K + 0.3*z;   % adiabate
    eta  = Ea/(R*T0_K);
    vALL = vALL .* exp( eta*( (T0_K./T_K) - 1) );
    % T-effect on viscosity; see Turcotte & Schubert, Geodynamics
    % Factor by which viscosity changes at every node
    vALL = min(vALL,1e8);
    vALL = max(vALL,1e-8);
    % limit viscosity change
end

% Pressure (depth) effect on viscosity (see Hirth & Kohlstedt 2003)
% =================================================================
vALL2 = vALL;
Va    = PHYSICS.Va; % m3/mol; activation volume
if Va>0
    P = abs(PHYSICS.g * Dens0 * NUMSCALE.L0 * depth(:));
        % lithostatic pressure at z
    vALL2 = vALL2 .* exp( P.*Va ./ (R.*T0_K) );
        % P-effect on viscosity (factor)
end

% Viscosity increase in lower mantle
% ==================================
vALL3 = vALL;
if isfield(PHYSICS,'a_UM2LM') && PHYSICS.a_UM2LM>0
    vALL3 = vALL3 * PHYSICS.a_UM2LM;
end

% Compute the minimum between pressure (depth) effect on viscosity and viscosity increase in the lower mantle 
% ===========================================================================================================
vALL = min([vALL2 vALL3],[],2); % This is done to avoid the viscosity "jump" in the upper mantle - lowe mantle transition

% New NODAL viscosity
% ===================
Visc = vALL .* Visc0 * NUMSCALE.Visc0;

% Viscosity cut-offs (if defined)
% ===============================
if isfield(PHYSICS,'maxVisc') && ~isempty(PHYSICS.maxVisc)
    Visc = min(Visc,PHYSICS.maxVisc*NUMSCALE.Visc0); % upper viscosity cut-off
end
if isfield(PHYSICS,'minVisc') && ~isempty(PHYSICS.minVisc)
    Visc = max(Visc,PHYSICS.minVisc*NUMSCALE.Visc0); % lower viscosity cut-off
end

ylim = 1e3 * [MESH.r_cmb MESH.r_surf] ./ NUMSCALE.L0;
semilogx(Visc(p),r(p),'b.-')
set(gca,'YLim',ylim);
set(gca,'XLim',[PHYSICS.minVisc*NUMSCALE.Visc0 PHYSICS.maxVisc*NUMSCALE.Visc0*10]);
xlabel('Visc (Pa-s)');
ylabel(sprintf('Distance (%s)',NUMSCALE.unit_L));
grid on

end

% #########################################################################

function plot_density(T,nodes,MESH,PHYSICS,NUMSCALE)

% Reference values
% =================
Dens0  = PHYSICS.Dens;
T0     = PHYSICS.T_ref;
r      = MESH.GCOORD_SPH(3,nodes);

% Temperature effect on density
% =============================
% (alpha is thermal expansion coeff, VAR.Dens is defined at the nodes)
Dens = Dens0 * (1 - PHYSICS.alpha*(T-T0)); % new density

% Density linearly varying with depth for the top half and constant for the bottom half 
% a = (Dens0 - Dens0*(1.1))/(6371-4921);
% b = Dens0 - a*6371;
% VAR.Dens = a*radius' + b;
% VAR.Dens(radius<4921) = 1.1*Dens0;

% Density cut-offs (if defined)
% =============================
if isfield(PHYSICS,'maxDens') && ~isempty(PHYSICS.maxDens)
    Dens = min(VAR.Dens,PHYSICS.maxDens); % upper density cut-off
end
if isfield(PHYSICS,'minDens') && ~isempty(PHYSICS.minDens)
    Dens = max(VAR.Dens,PHYSICS.minDens); % lower density cut-off
end

ylim = 1e3 * [MESH.r_cmb MESH.r_surf] ./ NUMSCALE.L0;
plot(Dens,r,'g.')
set(gca,'YLim',ylim);
xlabel('\rho (kg/m^3)');
ylabel(sprintf('Distance (%s)',NUMSCALE.unit_L));
grid on

end

% #########################################################################

function nodes_inside_region = nodes_inside_refined_region(MESH,lat_lon_region)

gTH_PT = MESH.GCOORD_SPH;
r_cmb  = MESH.r_cmb;
r_surf = MESH.r_surf;

theta_region               = (90 - lat_lon_region(1,:))*pi/180; % colatitude in radians
phi_region                 = lat_lon_region(2,:)*pi/180;        % longitude in radians
phi_region(phi_region < 0) = phi_region(phi_region < 0) + 2*pi; % values between 0 and 2pi
GCOORD_SPH_region_cmb      = [theta_region;phi_region;repmat(r_cmb,1,size(lat_lon_region,2))];
GCOORD_SPH_region_surf     = [theta_region;phi_region;repmat(r_surf,1,size(lat_lon_region,2))];
GCOORD_SPH_region          = [GCOORD_SPH_region_cmb GCOORD_SPH_region_surf];
EL2NOD_region              = delaunay(GCOORD_SPH_region');
EL2NOD_region              = uint32(EL2NOD_region');
[els,~,~]                  = tsearch2(GCOORD_SPH_region,EL2NOD_region,gTH_PT,[],[]);
nodes_inside_region        = els > 0;

% figure(344); clf
% GCOORD_SPH_plot      = GCOORD_SPH_region_rot;
% GCOORD_SPH_plot(3,:) = GCOORD_SPH_plot(3,:)/1000;
% hold on
% axis equal
% axis([0 pi 0 2*pi 3 6.5])
% view(142.5,30)
% grid on
% tetramesh(EL2NOD_region(1:4,:)',GCOORD_SPH_plot','FaceColor',[0.8 0 0],'FaceAlpha',0.3)
% scatter3(gTH_PT(1,:),gTH_PT(2,:),gTH_PT(3,:)/1000,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0.8])
% xlabel('\theta (rad)')
% ylabel('\phi (rad)')
% zlabel('r (10^3 km)')

end % END OF SUBFUNCTION nodes_inside_refined_region

% #########################################################################

function nodes_inside_craton = inside_craton(MESH,lat_lon_craton)

gTH_PT = MESH.GCOORD_SPH;
r_cmb  = MESH.r_cmb;
r_surf = MESH.r_surf;

theta_craton               = (90 - lat_lon_craton(1,:))*pi/180; % colatitude in radians
phi_craton                 = lat_lon_craton(2,:)*pi/180;        % longitude in radians
phi_craton(phi_craton < 0) = phi_craton(phi_craton < 0) + 2*pi; % values between 0 and 2pi
GCOORD_SPH_craton_cmb      = [theta_craton;phi_craton;repmat(r_cmb,1,size(lat_lon_craton,2))];
GCOORD_SPH_craton_surf     = [theta_craton;phi_craton;repmat(r_surf,1,size(lat_lon_craton,2))];
GCOORD_SPH_craton          = [GCOORD_SPH_craton_cmb GCOORD_SPH_craton_surf];

GCOORD_craton              = spherical2cartesian(GCOORD_SPH_craton);
GCOORD_craton_rot          = (MESH.RR2 * MESH.RR1)' * GCOORD_craton; % rotate coordinates from GPlates frame
GCOORD_SPH_craton_rot      = cartesian2spherical(GCOORD_craton_rot);
r_craton                   = [repmat(r_cmb,1,size(lat_lon_craton,2)) repmat(r_surf,1,size(lat_lon_craton,2))];
GCOORD_SPH_craton_rot(3,:) = r_craton;

EL2NOD_craton              = delaunay(GCOORD_SPH_craton_rot');
EL2NOD_craton              = uint32(EL2NOD_craton');
[els,~,~]                  = tsearch2(GCOORD_SPH_craton_rot,EL2NOD_craton,gTH_PT,[],[]);
nodes_inside_craton        = els > 0;

% figure(345); clf
% GCOORD_SPH_plot      = GCOORD_SPH_craton_rot;
% GCOORD_SPH_plot(3,:) = GCOORD_SPH_plot(3,:)/1000;
% hold on
% axis equal
% axis([0 pi 0 2*pi 3 6.5])
% view(142.5,30)
% grid on
% tetramesh(EL2NOD_craton(1:4,:)',GCOORD_SPH_plot','FaceColor',[0.8 0 0],'FaceAlpha',0.3)
% scatter3(gTH_PT(1,:),gTH_PT(2,:),gTH_PT(3,:)/1000,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0.8])
% xlabel('\theta (rad)')
% ylabel('\phi (rad)')
% zlabel('r (10^3 km)')

end % END OF SUBFUNCTION inside_craton