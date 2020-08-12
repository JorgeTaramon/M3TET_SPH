function [UBC,TBC,VAR] = bc_sph(MESH,COMM,SETTINGS,PHYSICS,VAR,NUMSCALE)
% Usage: [UBC,TBC,VAR] = bc_sph(MESH,COMM,SETTINGS,PHYSICS,VAR,NUMSCALE)
%
% Purpose: Define boundary conditions on temperature and velocity problem
%
% Input:
%   MESH      : [structure] : structure containing FE mesh data
%   COMM      : [structure] : inter-subdomain communication data
%   SETTINGS  : [structure] : model parameters
%   PHYSICS   : [structure] : structure containing physical properties
%   VAR       : [structure] : structure containing all major variables
%   NUMSCALE  : [structure] : numerical scaling parameters
%
% Output:
%   UBC       : [structure] : structure containing velocity boundary conditions
%   TBC       : [structure] : structure containing temperature boundary conditions
%   VAR       : [structure] : structure containing all major variables
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH March 2013
%

% All nodes at core-mantle boundary have PointID==301
% All nodes on top of the sphere have PointID==306

% Some variabels used below (inside the SUBFUNCTIONS):
%
% nUfix :: number of prescribed velocities degrees of freedom (dofs)
% iUfix :: indices of prescribed velocity dofs
% vUfix :: contains the velocities associated with iUfix
% nTfix :: number of prescribed temperatures
% iTfix :: indices of prescribed temperature nodes
% vTfix :: contains the temperatures associated iTfix
% nFlux :: number of prescribed heat flux boundary conditions
% iFlux :: indices of prescribed heat flux nodes
% vFlux :: contains the heat flux values associated with iFlux

switch SETTINGS.plate_model
    case 'WJM'
        [UBC,TBC,VAR] = plate_model_WJM(MESH,COMM,SETTINGS,PHYSICS,VAR);
        U             = zeros(3*MESH.nnod,1);
        U(UBC.ifix)   = UBC.vfix;
        VAR.Ur        = U(1:3:end);
        VAR.Ue        = U(2:3:end);
        VAR.Un        = U(3:3:end);
    case 'GPlates'
        [UBC,TBC,VAR] = plate_model_GPlates(MESH,COMM,PHYSICS,VAR);
        U             = zeros(3*MESH.nnod,1);
        U(UBC.ifix)   = UBC.vfix;
        VAR.Uth       = U(1:3:end);
        VAR.Uph       = U(2:3:end);
        VAR.Ur        = U(3:3:end);
    case 'Debug'
        [UBC,TBC,VAR] = plate_model_Debug(MESH,COMM,PHYSICS,VAR);
    otherwise
        error('SETTINGS.plate_model must be "WJM", "GPlates" or "Debug"')
end

return

% NEXT LINES ARE FOR TESTING MATRIX "MESH.RR_xyz2rEN":
omega = [0;0;1]; % test 1 (working)
% omega = [0;-1;0]; % test 2 (working)
% omega = [1;0;0]; % test 3
xyz   = MESH.GCOORD;
rmax  = max(sqrt(sum(xyz.^2)));
xyz   = xyz ./ rmax;
U_xyz = 100 * cross(repmat(omega,1,MESH.nnod),xyz);
U_xyz = U_xyz(:);

FigNo    = 60;
ind_plot = 6;
view_def = [-10 10];
% plot_domain_surfaces(FigNo  ,MESH,ind_plot,xyz(1,:));title('X');view(view_def);
% plot_domain_surfaces(FigNo+1,MESH,ind_plot,xyz(2,:));title('Y');view(view_def);
% plot_domain_surfaces(FigNo+2,MESH,ind_plot,xyz(3,:));title('Z');view(view_def);
plot_domain_surfaces(FigNo  ,MESH,ind_plot,U_xyz(1:3:end));title('Vel X');view(view_def);
plot_domain_surfaces(FigNo+1,MESH,ind_plot,U_xyz(2:3:end));title('Vel Y');view(view_def);
plot_domain_surfaces(FigNo+2,MESH,ind_plot,U_xyz(3:3:end));title('Vel Z');view(view_def);

U_rEN = MESH.RR_xyz2rEN * U_xyz;

plot_domain_surfaces(FigNo+3,MESH,ind_plot,U_rEN(1:3:end));title('Vel radial');view(view_def);
plot_domain_surfaces(FigNo+4,MESH,ind_plot,U_rEN(2:3:end));title('Vel East');view(view_def);
plot_domain_surfaces(FigNo+5,MESH,ind_plot,U_rEN(3:3:end));title('Vel North');view(view_def);

U_ThPhR = MESH.RR_new' * U_xyz;

plot_domain_surfaces(FigNo+6,MESH,ind_plot,U_ThPhR(1:3:end));title('Vel colatitudinal');view(view_def);
plot_domain_surfaces(FigNo+7,MESH,ind_plot,U_ThPhR(2:3:end));title('Vel longitudinal');view(view_def);
plot_domain_surfaces(FigNo+8,MESH,ind_plot,U_ThPhR(3:3:end));title('Vel radial');view(view_def);


U_colatitude = VAR.Un;
U_colatitude(DBnods(DBindx == 306)) = PHYSICS.Uth_Uph(1,:);
U_longitude = VAR.Ue;
U_longitude(DBnods(DBindx == 306))  = PHYSICS.Uth_Uph(2,:);
plot_domain_surfaces(FigNo+8,MESH,ind_plot,U_colatitude);title('Vel colatitudinal');view(view_def);
plot_domain_surfaces(FigNo+9,MESH,ind_plot,U_longitude);title('Vel longitudinal');view(view_def);


nods_cmb = find(MESH.PointID==301);
dofs_cmb = [3*nods_cmb-2 3*nods_cmb-1 3*nods_cmb]; % no slip core
% dofs_cmb = 3*nods_cmb-2; % free silp core
vals_cmb = zeros(1,length(dofs_cmb));
nods_top = find(MESH.PointID==306);
dofs_top = [3*nods_top-2 3*nods_top-1 3*nods_top];
vals_top = U_rEN(dofs_top)'; % spherical coordinates
% vals_top = U_xyz(dofs_top)'; % cartesian coordinates

UBC.ifix = [dofs_cmb dofs_top];
UBC.vfix = [vals_cmb vals_top];
UBC.nfix = length(UBC.ifix);

VAR.Ux   = U_xyz(1:3:end);
VAR.Uy   = U_xyz(2:3:end);
VAR.Uz   = U_xyz(3:3:end);

% U_xyz = MESH.RR_xyz2rEN' * U_rEN;
% 
% plot_domain_surfaces(FigNo+6,MESH,ind_plot,U_xyz(1:3:end));title('Vel X');view(view_def);
% plot_domain_surfaces(FigNo+7,MESH,ind_plot,U_xyz(2:3:end));title('Vel Y');view(view_def);
% plot_domain_surfaces(FigNo+8,MESH,ind_plot,U_xyz(3:3:end));title('Vel Z');view(view_def);

end % END OF FUNCTION bc_sph

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function [UBC,TBC,VAR] = plate_model_WJM(MESH,COMM,SETTINGS,PHYSICS,VAR)

% Make sure that all domain boundary indices, on which boundary conditions
% will be imposed, have been located in the mesh
indx    = [301 306];
DBnods  = prepare_DBnods(MESH,COMM,indx);
DBindx  = MESH.PointID(DBnods);
nDBnods = length(DBnods);

% Allocate arrays
iUbc    = 0; % counter for velocity boundary conditions
iUfix   = zeros(1,3*nDBnods); % list of prescribed velocity dofs
vUfix   = zeros(1,3*nDBnods); % list of the prescribed velocity values
iTbc    = 0; % counter for temperature boundary conditions
iTfix   = zeros(1,nDBnods); % list of prescribed temperature dofs (same as node number)
vTfix   = zeros(1,nDBnods); % list of the prescribed temperature values
iFbc    = 0; % counter for flux boundary conditions
iFlux   = zeros(1,nDBnods); % list of prescribed heat flux dofs (same as node number)
vFlux   = zeros(1,nDBnods); % list of the prescribed heat flux values

PlateInfo = zeros(nDBnods,3); itop = 0;

% Loop over all nodes on domain boundary
for ii=1:nDBnods
    node      = DBnods(ii); % global node number (= global temperature dof)
    ndf_r     = 3*node-2;   % global 1st velocity dof (r/vertical-direction)
    ndf_E     = 3*node-1;   % global 2nd velocity dof (E-direction)
    ndf_N     = 3*node;     % global 3rd velocity dof (N-direction)
    
    switch DBindx(ii)
        case 301 % CORE-MANTLE-BOUNDARY
            iTbc        = iTbc + 1;
            iTfix(iTbc) = node;
            vTfix(iTbc) = PHYSICS.Tbot;
            
            iUbc        = iUbc+1;
            iUfix(iUbc) = ndf_r; % this is the rotated
            vUfix(iUbc) = 0;
            
%             % --------------------------
%             % Uncomment these to make CMB be no-slip BC
%             % (comment these lines to make free slip core)
%             iUbc        = iUbc+1;
%             iUfix(iUbc) = ndf_N; % North-comp of Vel
%             vUfix(iUbc) = 0;
%             
%             iUbc        = iUbc+1;
%             iUfix(iUbc) = ndf_E; % East-comp of Vel
%             vUfix(iUbc) = 0;
%             % --------------------------
            
        case 306 % SURFACE
            iTbc        = iTbc + 1;
            iTfix(iTbc) = node;
            vTfix(iTbc) = PHYSICS.Ttop;
            
            % find out which plate the point is in
            % whichplate.m expects points to be unit vectors (?) stored as
            % x,y,z instead of lat,lon
            [platenum,~] = whichplateV(MESH.GCOORD(:,node),PHYSICS.PLATES);
            
            %             %                'jf' 'co' 'ca' 'ph'
            %             plate_delete  = [ 13   12   11   9  ];
            %             %                'pa' 'nz' 'na' 'pa'
            %             plate_instead = [  1    8    3    1 ];
            %             ind = find(platenum==plate_delete);
            %             if ~isempty(ind)
            %                 platenum = plate_instead(ind);
            %             end
            
            %             if platenum==5 % Antarctica
            %                 x = MESH.GCOORD(1,node);
            %                 y = MESH.GCOORD(2,node);
            %                 z = MESH.GCOORD(3,node);
            %
            % %                 figure(77);scatter3(x,y,z,10,'k','filled');hold on
            %                 % Find region around South pole
            %                 if z<-5800 && sqrt(x.^2+y.^2)<1500
            % %                     scatter3(x,y,z,10,'r','filled');hold on
            %                     continue
            %                 end
            %             end
            
            % calculate the cartesian velocity at this point (both as velEN and velXYZ)
            [velEN,~] = platevelocity(MESH.GCOORD(:,node),platenum,PHYSICS.PLATE_VELMODEL);
%             velEN(:) = 0;
            
            %             % CHECK:
            %             RR_nod = full(MESH.RR_xyz2rNE([ndf_r ndf_E ndf_N],[ndf_r ndf_E ndf_N]));
            %             check  = [0;velEN] - RR_nod * velXYZ;
            %             if max(abs(check>1e-12))
            %                 error('Rotation matrix corrupt?');
            %             end
            
            %             VAR.Ux(node) = velXYZ(1);
            %             VAR.Uy(node) = velXYZ(2);
            %             VAR.Uz(node) = velXYZ(3);
            %             VAR.Ur(node) = 0;
            %             VAR.Ue(node) = velEN(1);
            %             VAR.Un(node) = velEN(2);
            
            % find out which Nataf Region this point is in (cratonic
            % lithosphere if NatafNum = 1 or 2, young cont = 3, oceans =
            % 4,5)
            [NatafNum] = whichNatafRegion(MESH.GCOORD(:,node),PHYSICS.NATAF_CRATON);
            
            itop              = itop + 1;
            PlateInfo(itop,1) = node;
            PlateInfo(itop,2) = platenum;
            PlateInfo(itop,3) = NatafNum; % store Nataf Tect Regionalization number for this point -- so can be used later for plotting
            
            set_plate_motion = 0;
            if NatafNum <= 2 || SETTINGS.useNataf==0
                set_plate_motion = 1;
            end
            
            if set_plate_motion
                % Set plate motion defined in r-, E-, N-coordinates
                iUbc        = iUbc+1;
                iUfix(iUbc) = ndf_E; % East-comp of Vel
                vUfix(iUbc) = velEN(1);
                %                 vUfix(iUbc) = velXYZ(2);
                
                iUbc        = iUbc+1;
                iUfix(iUbc) = ndf_N; % North-comp of Vel
                vUfix(iUbc) = velEN(2);
                %                 vUfix(iUbc) = velXYZ(1);
                
                iUbc        = iUbc+1;
                iUfix(iUbc) = ndf_r; % radial-comp of Vel
                vUfix(iUbc) = 0;
                %                 vUfix(iUbc) = velXYZ(3);
                
            else
                % For FREE-SLIP node on surface (ALSO NEED TO INCLUDE IN ROTATED local coords
                iUbc        = iUbc+1;
                iUfix(iUbc) = ndf_r; % this is the rotated V_r dof
                vUfix(iUbc) = 0;
                
                %                 % possible boundary condition is the plate speed
                %                 % fixing everything to zero for now
                %                 iUbc        = iUbc+1;
                %                 iUfix(iUbc) = ndf_x; % x-comp of Vel
                %      %           vUfix(iUbc) = velXYZ(1);
                %                 vUfix(iUbc) = 0.; % debug
                %
                %                 iUbc        = iUbc+1;
                %                 iUfix(iUbc) = ndf_y; % y-comp of Vel
                %      %           vUfix(iUbc) = velXYZ(2);
                %                 vUfix(iUbc) = 0.;
                %
                %                 iUbc        = iUbc+1;
                %                 iUfix(iUbc) = ndf_z; % z-comp of Vel
                %      %           vUfix(iUbc) = velXYZ(3);
                %                 vUfix(iUbc) = 0.;
            end
            
        otherwise
            error(' Unknown boundary index for a sphere.');
    end
end

% shrink arrays to correct size
iUfix(iUbc+1:end) = [];
vUfix(iUbc+1:end) = [];
iTfix(iTbc+1:end) = [];
vTfix(iTbc+1:end) = [];
iFlux(iFbc+1:end) = [];
vFlux(iFbc+1:end) = [];
PlateInfo(itop+1:end,:) = [];

% Initialize boundary condition structures
UBC = struct('nfix' , iUbc,...
             'ifix' , iUfix,...
             'vfix' , vUfix,...
             'PlateInfo',PlateInfo);
TBC = struct('nfix' , iTbc,...
             'ifix' , iTfix,...
             'vfix' , vTfix,...
             'nflux', iFbc,...
             'iflux', iFlux,...
             'vflux', vFlux );

VAR.T(iTfix) = vTfix; % update temperature field to include BCs

VAR.Plate                     = zeros(MESH.nnod,1);
VAR.Plate(UBC.PlateInfo(:,1)) = UBC.PlateInfo(:,2);

end % END OF SUBFUNCTION plate_model_WJM
        
% #########################################################################     

function [UBC,TBC,VAR] = plate_model_GPlates(MESH,COMM,PHYSICS,VAR)

% load('MESH_GPLATES_BCs.mat');
% u_sph_rot1 = interpolate_GPlates_BCs_method1(MESH,MESH_GPLATES_BCs,PHYSICS,COMM); % WORSE RESULTS 
% u_sph_rot2 = interpolate_GPlates_BCs_method2(MESH,MESH_GPLATES_BCs,PHYSICS,COMM); % BETTER RESUTLS 

% COMPUTE TEMPERATURE AND VELOCITY BCS FOR CMB AND SURFACE
% ========================================================
% Make sure that all domain boundary indices, on which boundary conditions will be imposed, have been located in the mesh
indx    = [301 306];
DBnods  = prepare_DBnods(MESH,COMM,indx);
DBindx  = MESH.PointID(DBnods);
nDBnods = length(DBnods);

% Allocate arrays
iUbc      = 0;                  % counter for velocity boundary conditions
iUfix     = zeros(1,3*nDBnods); % list of prescribed velocity dofs
vUfix     = zeros(1,3*nDBnods); % list of the prescribed velocity values
iTbc      = 0;                  % counter for temperature boundary conditions
iTfix     = zeros(1,nDBnods);   % list of prescribed temperature dofs (same as node number)
vTfix     = zeros(1,nDBnods);   % list of the prescribed temperature values
iFbc      = 0;                  % counter for flux boundary conditions
iFlux     = zeros(1,nDBnods);   % list of prescribed heat flux dofs (same as node number)
vFlux     = zeros(1,nDBnods);   % list of the prescribed heat flux values
iUth      = 0;                  % counter for reading Uth_Uph matrix
iUph      = 0;                  % counter for reading Uth_Uph matrix

% Rotate GPlates velocities (actual Earth) into the spherical mesh frame (where the refined region is centered at theta = 90, phi = 90)
sph             = PHYSICS.V_BCs(1).latlon;
sph(1,:)        = 90 - sph(1,:); % transform latitude to colatitude
sph             = (pi/180)*sph;  % transform degrees to radians
sph             = [sph; MESH.r_surf*ones(1,size(sph,2))];
cart            = spherical2cartesian(sph);
TT              = make_rotation_matrix(cart); % transformation matrix between spherical and Cartesian coordinates
u_sph           = [PHYSICS.V_BCs(1).Uth_Uph; zeros(1,size(PHYSICS.V_BCs(1).Uth_Uph,2))];
u_cart          = TT * u_sph(:);
u_cart_rot      = (MESH.RR2 * MESH.RR1)' * reshape(u_cart,3,[]); % rotate velocities
cart_rot        = (MESH.RR2 * MESH.RR1)' * cart;  % rotate coordinates
TT_rot          = make_rotation_matrix(cart_rot); % transformation matrix between spherical and Cartesian rotated coordinates
u_sph_rot       = TT_rot'* u_cart_rot(:);
u_sph_rot       = reshape(u_sph_rot,3,[]);
u_sph_rot(3,:)  = [];
PHYSICS.V_BCs(1).Uth_Uph = u_sph_rot;

% Loop over all nodes on domain boundary
for ii=1:nDBnods
    node   = DBnods(ii); % global node number (= global temperature dof)
    ndf_th = 3*node-2;   % global 1st velocity dof (colatitude (theta) direction)
    ndf_ph = 3*node-1;   % global 2nd velocity dof (longitude (phi) direction)
    ndf_r  = 3*node;     % global 3rd velocity dof (radial (r) direction)
    
    switch DBindx(ii)
        case 301 % CORE-MANTLE-BOUNDARY
            % Temperature BC
            iTbc        = iTbc + 1;
            iTfix(iTbc) = node;
            vTfix(iTbc) = PHYSICS.Tbot;
            
            % Velocity BC (free slip)
            iUbc        = iUbc+1;
            iUfix(iUbc) = ndf_r; % this is the rotated
            vUfix(iUbc) = 0;
            
%             % --------------------------
%             % Uncomment these to make CMB be no-slip BC
%             % (comment these lines to make free slip core)
%             iUbc        = iUbc+1;
%             iUfix(iUbc) = ndf_th; % colatitude-comp of Vel (North-comp)
%             vUfix(iUbc) = 0;
%             
%             iUbc        = iUbc+1;
%             iUfix(iUbc) = ndf_ph; % longitude-comp of Vel (East-comp)
%             vUfix(iUbc) = 0;
%             % --------------------------
            
        case 306 % SURFACE
            % Temperature BC
            iTbc        = iTbc + 1;
            iTfix(iTbc) = node;
            vTfix(iTbc) = PHYSICS.Ttop;
            
            % Velocity BC
            iUth        = iUth+1;
            iUbc        = iUbc+1;
            iUfix(iUbc) = ndf_th; % dof colatitude-comp of Vel (number)
            vUfix(iUbc) = PHYSICS.V_BCs(1).Uth_Uph(1,iUth); % colatitude-comp of Vel (value)
            
            iUph        = iUph+1;
            iUbc        = iUbc+1;
            iUfix(iUbc) = ndf_ph; % dof longitude-comp of Vel (number)
            vUfix(iUbc) = PHYSICS.V_BCs(1).Uth_Uph(2,iUph); % longitude-comp of Vel (value)
            
            iUbc        = iUbc+1;
            iUfix(iUbc) = ndf_r; % this is the rotated
            vUfix(iUbc) = 0;     % tangent motion to the surface
    end
end

% shrink arrays to correct size
iUfix(iUbc+1:end)       = [];
vUfix(iUbc+1:end)       = [];
iTfix(iTbc+1:end)       = [];
vTfix(iTbc+1:end)       = [];
iFlux(iFbc+1:end)       = [];
vFlux(iFbc+1:end)       = [];

if ~isempty(PHYSICS.PLUME)
    % COMPUTE TEMPERATURE AND VELOCITY BCS FOR A PLUME
    % ================================================
    T_plume                 = PHYSICS.PLUME.T_plume;
    colat_lon_plume_center  = PHYSICS.PLUME.colat_lon_plume_center; % center of the plume (colatitude,longitude) [degrees]
    r_plume                 = PHYSICS.PLUME.r_plume; % radius of the plume [km]
    l_TBC                   = PHYSICS.PLUME.l_TBC;   % length of the plume [km] where Plume Temp BCs are imposed
    d_TBC                   = PHYSICS.PLUME.d_TBC;   % plume bottom depth [km] for imposing Temp Bcs
    l_UBC                   = PHYSICS.PLUME.l_UBC;   % length of the plume [km] where Plume Vel BCs are imposed
    d_UBC                   = PHYSICS.PLUME.d_UBC;   % plume bottom depth [km] for imposing Vel Bcs
    u_max                   = PHYSICS.PLUME.u_max;    
    
    % nodes and Temperature for imposing plume TBC
    [nodes_inside_plume_TBC,~,T_nodes_inside_plume] = ...
        inside_plume(MESH,VAR.T,T_plume,colat_lon_plume_center,r_plume,l_TBC,d_TBC,u_max);
    DBnods_plume_TBC  = find(nodes_inside_plume_TBC);
    nDBnods_plume_TBC = length(DBnods_plume_TBC);
    
    % nodes and velocities for imposing plume UBC
    [nodes_inside_plume_UBC,Ur_plume,~] = ...
        inside_plume(MESH,VAR.T,T_plume,colat_lon_plume_center,r_plume,l_UBC,d_UBC,u_max);
    DBnods_plume_UBC = find(nodes_inside_plume_UBC);
    
    % Allocate arrays
    iUbc_plume      = 0;                            % counter for velocity boundary conditions
    iUfix_plume     = zeros(1,3*nDBnods_plume_TBC); % list of prescribed velocity dofs
    vUfix_plume     = zeros(1,3*nDBnods_plume_TBC); % list of the prescribed velocity values
    iTbc_plume      = 0;                            % counter for temperature boundary conditions
    iTfix_plume     = zeros(1,nDBnods_plume_TBC);   % list of prescribed temperature dofs (same as node number)
    vTfix_plume     = zeros(1,nDBnods_plume_TBC);   % list of the prescribed temperature values
    iUr             = 0;                            % counter for reading Ur_plume vector
    iT              = 0;                            % counter for reading T vector
    
    % Loop over all nodes inside the plume
    for i=1:nDBnods_plume_TBC
        node_plume   = DBnods_plume_TBC(i); % global node number (= global temperature dof)
        
        if ismember(node_plume,DBnods_plume_UBC) % node inside the top disc of the plume
            ndf_th_plume = 3*node_plume-2;  % global 1st velocity dof (colatitude (theta) direction)
            ndf_ph_plume = 3*node_plume-1;  % global 2nd velocity dof (longitude (phi) direction)
            ndf_r_plume  = 3*node_plume;    % global 3rd velocity dof (radial (r) direction)
            
            % Temperature BC
            iT                      = iT+1;
            iTbc_plume              = iTbc_plume + 1;
            iTfix_plume(iTbc_plume) = node_plume;
            vTfix_plume(iTbc_plume) = T_nodes_inside_plume(iT);
            
            % Velocity BC
%             iUbc_plume              = iUbc_plume+1;
%             iUfix_plume(iUbc_plume) = ndf_th_plume;  % dof colatitude-comp of Vel (number)
%             vUfix_plume(iUbc_plume) = 0;             % colatitude-comp of Vel (value)
%             
%             iUbc_plume              = iUbc_plume+1;
%             iUfix_plume(iUbc_plume) = ndf_ph_plume;  % dof longitude-comp of Vel (number)
%             vUfix_plume(iUbc_plume) = 0;             % longitude-comp of Vel (value)
            
            iUr                     = iUr+1;
            iUbc_plume              = iUbc_plume+1;
            iUfix_plume(iUbc_plume) = ndf_r_plume;   % dof radial-comp of Vel (number)
            vUfix_plume(iUbc_plume) = Ur_plume(iUr); % radial-comp of Vel (value)
        else
            % Temperature BC
            iT                      = iT+1;
            iTbc_plume              = iTbc_plume + 1;
            iTfix_plume(iTbc_plume) = node_plume;
            vTfix_plume(iTbc_plume) = T_nodes_inside_plume(iT);
        end
    end
    
    % shrink arrays to correct size
    iUfix_plume(iUbc_plume+1:end)       = [];
    vUfix_plume(iUbc_plume+1:end)       = [];
    iTfix_plume(iTbc_plume+1:end)       = [];
    vTfix_plume(iTbc_plume+1:end)       = [];
    
    % Save plume velocity BCs to load them again when updating velocity BCs on surface for each time step 
    PLUME.iUbc_plume  = iUbc_plume;
    PLUME.iUfix_plume = iUfix_plume;
    PLUME.vUfix_plume = vUfix_plume;
else
    iUbc_plume  = 0;
    iUfix_plume = [];
    vUfix_plume = [];
    iTbc_plume  = 0;
    iTfix_plume = [];
    vTfix_plume = [];
    PLUME       = [];
end

% Initialize boundary condition structures
UBC = struct('nfix'      , iUbc + iUbc_plume,        ...
             'ifix'      , [iUfix iUfix_plume],      ...
             'vfix'      , [vUfix vUfix_plume],      ...
             'PlateInfo' , PHYSICS.V_BCs(1).PlateID, ...
             'PLUME'     , PLUME);
TBC = struct('nfix'      , iTbc + iTbc_plume,   ...
             'ifix'      , [iTfix iTfix_plume], ...
             'vfix'      , [vTfix vTfix_plume], ...
             'nflux'     , iFbc,                ...
             'iflux'     , iFlux,               ...
             'vflux'     , vFlux );

VAR.T([iTfix iTfix_plume])     = [vTfix vTfix_plume]; % update temperature field to include BCs
VAR.Plate                      = zeros(MESH.nnod,1);
VAR.Plate(MESH.PointID == 306) = UBC.PlateInfo;

end % END OF SUBFUNCTION plate_model_GPlates

% #########################################################################     

function [UBC,TBC,VAR] = plate_model_Debug(MESH,COMM,PHYSICS,VAR)

inod_cmb    = find(MESH.PointID==301);
inod_surf   = find(MESH.PointID==306);
% TEMPERATUE BOUNDARY CONDITION
TBC.ifix    = inod_cmb;
TBC.vfix    = PHYSICS.Tbot * ones(1,length(inod_cmb));
TBC.ifix    = [TBC.ifix inod_surf];
TBC.vfix    = [TBC.vfix PHYSICS.Ttop * ones(1,length(inod_surf))];
% update temperature field to include BCs
VAR.T(TBC.ifix) = TBC.vfix; 


% VELOCITY BOUNDARY CONDITION
idof_th     = 3*[inod_cmb inod_surf]-2;
idof_ph     = 3*[inod_cmb inod_surf]-1;
idof_r      = 3*[inod_cmb inod_surf];

%==========================================================================
% Zero velocity BC
UBC.ifix    = [idof_th idof_ph idof_r];
UBC.vfix    = [zeros(1,length(idof_th)) zeros(1,length(idof_ph)) zeros(1,length(idof_r))];

% %==========================================================================
% % Free slip (velocity is set zero at radial direction)
% UBC.ifix    = idof_r;
% UBC.vfix    = (zeros(1,length(idof_r)));

% %==========================================================================
% % Free slip (velocity is set zero at radial direction) except for two nodes (where velocity is set zero at all directions): 
% % - 1st fixed node at North Pole on the top surface (theta = 0°)
% % - 2nd fixed node at intersection of equator and phi = 180° meridian on the top surface (theta = 90°, phi = 180°)
% % This is to avoid a spinning flow around the X and Z axes
% GCOORD_SPH      = cartesian2spherical(MESH.GCOORD);
% GCOORD_SPH(1,:) = GCOORD_SPH(1,:)*180/pi; % transform theta to degrees
% GCOORD_SPH(2,:) = GCOORD_SPH(2,:)*180/pi; % transform phi to degrees
% node1           = inod_surf(GCOORD_SPH(1,inod_surf) < 5);
% idof_theta1     = 3*node1(1) - 2;
% idof_phi1       = 3*node1(1) - 1;
% vnod_theta1     = 0;
% vnod_phi1       = 0;
% node2           = inod_surf(GCOORD_SPH(1,inod_surf) >  85 & GCOORD_SPH(1,inod_surf) <  95 & ...
%                             GCOORD_SPH(2,inod_surf) > 185 & GCOORD_SPH(2,inod_surf) < 195);
% % idof_theta2     = 3*node2(1) - 2;
% idof_phi2       = 3*node2(1) - 1;
% % vnod_theta2     = 0;
% vnod_phi2       = 0;
% UBC.ifix        = [idof_theta1 idof_phi1 idof_phi2                idof_r  ];
% UBC.vfix        = [vnod_theta1 vnod_phi1 vnod_phi2 zeros(1,length(idof_r))];

% %==========================================================================
% % Zero velocity BC at surface and Free slip at CMB
% idof_th_cmb  = 3*inod_cmb-2;
% idof_ph_cmb  = 3*inod_cmb-1;
% idof_r_cmb   = 3*inod_cmb;
% idof_th_surf = 3*inod_surf-2;
% idof_ph_surf = 3*inod_surf-1;
% idof_r_surf  = 3*inod_surf;
% UBC.ifix    = [               idof_r_cmb                  idof_th_surf                  idof_ph_surf                  idof_r_surf];
% UBC.vfix    = [zeros(1,length(idof_r_cmb)) zeros(1,length(idof_th_surf)) zeros(1,length(idof_ph_surf)) zeros(1,length(idof_r_surf))];

% %==========================================================================
% % Prescribed velocities in theta direction (theta is positive in clockwise
% % direcction) for 2 plates:
% %   1st plate: theta = 60° - 120° and phi =  60° - 120° with velocity v_theta =  40 mm/yr (v_colatitudinal) 
% %   2nd plate: theta = 60° - 120° and phi = 240° - 300° with velocity v_theta =  40 mm/yr (v_colatitudinal)
% GCOORD_SPH      = cartesian2spherical(MESH.GCOORD);
% GCOORD_SPH(1,:) = GCOORD_SPH(1,:)*180/pi; % transform theta to degrees
% GCOORD_SPH(2,:) = GCOORD_SPH(2,:)*180/pi; % transform phi to degrees
% idof_theta1     = 3*(inod_surf(GCOORD_SPH(1,inod_surf) >  60 & GCOORD_SPH(1,inod_surf) < 120 & ...
%                                GCOORD_SPH(2,inod_surf) >  60 & GCOORD_SPH(2,inod_surf) < 120))-2;
% vnod_theta1     = 40*ones(1,length(idof_theta1));
% idof_phi1       = 3*(inod_surf(GCOORD_SPH(1,inod_surf) >  60 & GCOORD_SPH(1,inod_surf) < 120 & ...
%                                GCOORD_SPH(2,inod_surf) >  60 & GCOORD_SPH(2,inod_surf) < 120))-1;
% vnod_phi1       = zeros(1,length(idof_phi1));
% idof_theta2     = 3*(inod_surf(GCOORD_SPH(1,inod_surf) >  60 & GCOORD_SPH(1,inod_surf) < 120 & ...
%                                GCOORD_SPH(2,inod_surf) > 240 & GCOORD_SPH(2,inod_surf) < 300))-2;
% vnod_theta2     = 40*ones(1,length(idof_theta2));
% idof_phi2       = 3*(inod_surf(GCOORD_SPH(1,inod_surf) >  60 & GCOORD_SPH(1,inod_surf) < 120 & ...
%                                GCOORD_SPH(2,inod_surf) > 240 & GCOORD_SPH(2,inod_surf) < 300))-1;
% vnod_phi2       = zeros(1,length(idof_phi2));
% UBC.ifix        = [idof_theta1 idof_phi1 idof_theta2 idof_phi2                idof_r  ];
% UBC.vfix        = [vnod_theta1 vnod_phi1 vnod_theta2 vnod_phi2 zeros(1,length(idof_r))];
% FigNo    = 71; 
% close(figure(FigNo))
% ind_plot = 2;
% view_def = [142.5 30];
% U_colatitude = VAR.Uth;
% inod_surf1 = inod_surf(GCOORD_SPH(1,inod_surf) >  60 & GCOORD_SPH(1,inod_surf) < 120 & ...
%                        GCOORD_SPH(2,inod_surf) >  60 & GCOORD_SPH(2,inod_surf) < 120);
% inod_surf2 = inod_surf(GCOORD_SPH(1,inod_surf) >  60 & GCOORD_SPH(1,inod_surf) < 120 & ...
%                        GCOORD_SPH(2,inod_surf) > 240 & GCOORD_SPH(2,inod_surf) < 300);
% U_colatitude([inod_surf1 inod_surf2]) = [vnod_theta1 vnod_theta2];
% plot_domain_surfaces(FigNo,MESH,ind_plot,U_colatitude);title('Vel colatitudinal');view(view_def);

% %==========================================================================
% % Prescribed velocities for a rigid rotation using Z axis as Euler Pole
% r          = 6371;     % [km]     -> radius where vel is maximum
% v          = 80;       % [km/Myr] -> lineal speed at r
% omega      = v/r;      % rad/Myr
% E          = [0 0 1]'; % unit vector for Euler Pole (Z-axis)
% U_cart     = omega * cross(repmat(E,1,MESH.nnod),MESH.GCOORD);
% U_sph      = MESH.RR'*U_cart(:);
% U_sph      = reshape(U_sph,3,[]);
% U_sph      = U_sph(:,inod_surf);
% % Make sure that all domain boundary indices, on which boundary conditions will be imposed, have been located in the mesh
% indx    = [301 306];
% DBnods  = prepare_DBnods(MESH,COMM,indx);
% DBindx  = MESH.PointID(DBnods);
% nDBnods = length(DBnods);
% % Allocate arrays
% iUbc      = 0;                  % counter for velocity boundary conditions
% iUfix     = zeros(1,3*nDBnods); % list of prescribed velocity dofs
% vUfix     = zeros(1,3*nDBnods); % list of the prescribed velocity values
% iUth      = 0;                  % counter for reading Uth_Uph matrix
% iUph      = 0;                  % counter for reading Uth_Uph matrix
% % Loop over all nodes on domain boundary
% for ii=1:nDBnods
%     node   = DBnods(ii); % global node number (= global temperature dof)
%     ndf_th = 3*node-2;   % global 1st velocity dof (colatitude (theta) direction)
%     ndf_ph = 3*node-1;   % global 2nd velocity dof (longitude (phi) direction)
%     ndf_r  = 3*node;     % global 3rd velocity dof (radial (r) direction)
%     switch DBindx(ii)
%         case 301 % CORE-MANTLE-BOUNDARY
%             % Velocity BC (free slip)
%             iUbc        = iUbc+1;
%             iUfix(iUbc) = ndf_r; % this is the rotated
%             vUfix(iUbc) = 0;
%         case 306 % SURFACE
%             % Velocity BC (rotation around Z axis)
%             iUth        = iUth+1;
%             iUbc        = iUbc+1;
%             iUfix(iUbc) = ndf_th; % dof colatitude-comp of Vel (number)
%             vUfix(iUbc) = 0;      % colatitude-comp of Vel (value)
%             
%             iUph        = iUph+1;
%             iUbc        = iUbc+1;
%             iUfix(iUbc) = ndf_ph; % dof longitude-comp of Vel (number)
%             vUfix(iUbc) = U_sph(2,iUph); % longitude-comp of Vel (value)
%             
%             iUbc        = iUbc+1;
%             iUfix(iUbc) = ndf_r;
%             vUfix(iUbc) = 0;      % tangent motion to the surface
%     end
% end
% % shrink arrays to correct size
% iUfix(iUbc+1:end)       = [];
% vUfix(iUbc+1:end)       = [];
% 
% UBC.nfix = iUbc;
% UBC.ifix = iUfix;
% UBC.vfix = vUfix;
% %
% % FigNo    = 71; 
% % close(figure(FigNo))
% % ind_plot = 2;
% % view_def = [142.5 30];
% % speed    = sqrt(sum(U_cart.^2));
% % plot_domain_surfaces(FigNo,MESH,ind_plot,speed);title('Speed rigid rotation');view(view_def);
% % hold on
% % quiver3(MESH.GCOORD(1,:),MESH.GCOORD(2,:),MESH.GCOORD(3,:),U_cart(1,:),U_cart(2,:),U_cart(3,:),0)

% %==========================================================================
% % Prescribed velocities for 2 plates making a ridge and a trench
% r          = 6371;     % [km]     -> radius where vel is maximum
% v          = 80;       % [km/Myr] -> lineal speed at r
% E          = [0 0 1]'; % unit vector for Euler Pole (Z-axis)
% omega_CW   = v/r;      % rad/Myr
% U_cart_CW  = omega_CW * cross(repmat(E,1,MESH.nnod),MESH.GCOORD);
% U_sph_CW   = MESH.RR'*U_cart_CW(:);
% U_sph_CW   = reshape(U_sph_CW,3,[]);
% U_sph_CW   = U_sph_CW(:,inod_surf);
% omega_CCW  = -v/r;      % rad/Myr
% U_cart_CCW = omega_CCW * cross(repmat(E,1,MESH.nnod),MESH.GCOORD);
% U_sph_CCW  = MESH.RR'*U_cart_CCW(:);
% U_sph_CCW  = reshape(U_sph_CCW,3,[]);
% U_sph_CCW  = U_sph_CCW(:,inod_surf);
% % Make sure that all domain boundary indices, on which boundary conditions will be imposed, have been located in the mesh
% indx    = [301 306];
% DBnods  = prepare_DBnods(MESH,COMM,indx);
% DBindx  = MESH.PointID(DBnods);
% nDBnods = length(DBnods);
% % Allocate arrays
% iUbc      = 0;                  % counter for velocity boundary conditions
% iUfix     = zeros(1,3*nDBnods); % list of prescribed velocity dofs
% vUfix     = zeros(1,3*nDBnods); % list of the prescribed velocity values
% iUth      = 0;                  % counter for reading Uth_Uph matrix
% iUph      = 0;                  % counter for reading Uth_Uph matrix
% % Loop over all nodes on domain boundary
% for ii=1:nDBnods
%     node   = DBnods(ii); % global node number (= global temperature dof)
%     ndf_th = 3*node-2;   % global 1st velocity dof (colatitude (theta) direction)
%     ndf_ph = 3*node-1;   % global 2nd velocity dof (longitude (phi) direction)
%     ndf_r  = 3*node;     % global 3rd velocity dof (radial (r) direction)
%     switch DBindx(ii)
%         case 301 % CORE-MANTLE-BOUNDARY
%             % Velocity BC (free slip)
%             iUbc        = iUbc+1;
%             iUfix(iUbc) = ndf_r; % this is the rotated
%             vUfix(iUbc) = 0;
%         case 306 % SURFACE
%             % Velocity BC (rotation around Z axis)
%             iUth        = iUth+1;
%             iUbc        = iUbc+1;
%             iUfix(iUbc) = ndf_th; % dof colatitude-comp of Vel (number)
%             vUfix(iUbc) = 0;      % colatitude-comp of Vel (value)
%             
%             iUph        = iUph+1;
%             iUbc        = iUbc+1;
%             iUfix(iUbc) = ndf_ph; % dof longitude-comp of Vel (number)
%             if MESH.GCOORD_SPH(2,DBnods(ii)) >= pi/2 && MESH.GCOORD_SPH(2,DBnods(ii)) < 3*pi/2
%                 vUfix(iUbc) = U_sph_CW(2,iUph); % longitude-comp of Vel (value)
%             else
%                 vUfix(iUbc) = U_sph_CCW(2,iUph); % longitude-comp of Vel (value)
%             end
%             iUbc        = iUbc+1;
%             iUfix(iUbc) = ndf_r;
%             vUfix(iUbc) = 0;      % tangent motion to the surface
%     end
% end
% % shrink arrays to correct size
% iUfix(iUbc+1:end)       = [];
% vUfix(iUbc+1:end)       = [];
% 
% UBC.nfix = iUbc;
% UBC.ifix = iUfix;
% UBC.vfix = vUfix;

end % END OF SUBFUNCTION plate_model_Debug

% #########################################################################     

function [DBnods,indx2bc] = prepare_DBnods(MESH,COMM,indx)

% Make sure that all domain boundary indices, on which boundary conditions
% will be imposed, have been located in the mesh
indx      = indx';
DBnods    = find(MESH.PointID);
DBindices = unique(MESH.PointID(DBnods)); % make a unique list of domain boundary indices
if COMM.nsd>1
    DBindices_par = COMM.unique(DBindices);
    indx_par      = COMM.unique(indx);
    missing       = setdiff(indx_par,DBindices_par);
else
    missing = setdiff(indx,DBindices);
end
if ~isempty(missing)
    error(['The following indices are in setbc but have no associated index in DBnods: '...
           num2str(missing(:)')])
end
if COMM.nsd>1
    missing  = setdiff(DBindices_par,indx_par);
else
    missing = setdiff(DBindices,indx);
end
if ~isempty(missing)
    disp(['WARNING: Boundaries with the following indices have no boundary condition: '...
           num2str(missing(:)')])
end

% Pointer from domain boundary index to the corresponding row in setbc 
% and bcval
indx2bc(indx) = 1:length(indx);

end % END OF SUBFUNCTION prepare_DBnods

% #########################################################################

function [nodes_inside_plume,u_plume,T_nodes_inside_plume] = inside_plume(MESH,T,T_plume,colat_lon_plume_center,r_plume,l_plume,d_plume,u_max)

gTH_PT         = MESH.GCOORD_SPH;
r_plume_bottom = MESH.r_surf - d_plume;    % radial distance of the bottom of the plume
r_plume_top    = r_plume_bottom + l_plume; % radial distance of the top of the plume

theta_plume_center                     = colat_lon_plume_center(1)*pi/180; % colatitude in radians
phi_plume_center                       = colat_lon_plume_center(2)*pi/180; % longitude in radians
phi_plume_center(phi_plume_center < 0) = phi_plume_center(phi_plume_center < 0) + 2*pi; % values between 0 and 2pi

GCOORD_plume_center         = spherical2cartesian([theta_plume_center; phi_plume_center; MESH.r_surf]);
GCOORD_plume_center_rot     = (MESH.RR2 * MESH.RR1)' * GCOORD_plume_center; % rotate coordinates of the plume center from GPlates frame
GCOORD_SPH_plume_center_rot = cartesian2spherical(GCOORD_plume_center_rot);

% % % convert radius of the plume (in km) into degrees and then into radians at the depth r_plume_bottom 
% % r_plume_disc_deg = r_plume / (pi*r_plume_bottom/180);
% % r_plume_disc_rad = r_plume_disc_deg*pi/180;
% % 
% % % create points in the tail boundary
% % alpha       = (0:5:355)*pi/180;
% % theta_plume = zeros(1,size(alpha,2));
% % phi_plume   = zeros(1,size(alpha,2));
% % for i = 1:size(alpha,2)
% %    theta_plume(i) = GCOORD_SPH_plume_center_rot(1) + r_plume_disc_rad * cos(alpha(i));
% %    phi_plume(i)   = GCOORD_SPH_plume_center_rot(2) + r_plume_disc_rad * sin(alpha(i));
% % end
% % 
% % GCOORD_SPH_plume_bottom = [theta_plume;phi_plume;repmat(r_plume_bottom,1,size(alpha,2))];
% % GCOORD_SPH_plume_top    = [theta_plume;phi_plume;repmat(r_plume_top,1,size(alpha,2))];
% % GCOORD_SPH_plume        = [GCOORD_SPH_plume_bottom GCOORD_SPH_plume_top];
% % 
% % EL2NOD_plume            = delaunay(GCOORD_SPH_plume');
% % EL2NOD_plume            = uint32(EL2NOD_plume');
% % [els,~,~]               = tsearch2(GCOORD_SPH_plume,EL2NOD_plume,gTH_PT,[],[]);
% % nodes_inside_plume      = els > 0;
% % 
% % figure(347); clf
% % GCOORD_SPH_plot      = GCOORD_SPH_plume;
% % GCOORD_SPH_plot(3,:) = GCOORD_SPH_plot(3,:)/1000;
% % hold on
% % axis equal
% % % axis([0 pi 0 2*pi 3 6.5])
% % view(142.5,30)
% % grid on
% % tetramesh(EL2NOD_plume(1:4,:)',GCOORD_SPH_plot','FaceColor',[0.8 0 0],'FaceAlpha',0.3)
% % scatter3(MESH.GCOORD_SPH(1,nodes_inside_plume), ...
% %          MESH.GCOORD_SPH(2,nodes_inside_plume), ...
% %          MESH.GCOORD_SPH(3,nodes_inside_plume)/1000,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0.8])
% % xlabel('\theta (rad)')
% % ylabel('\phi (rad)')
% % zlabel('r (10^3 km)')
% % 
% % theta_nodes_inside_plume = gTH_PT(1,nodes_inside_plume);
% % phi_nodes_inside_plume   = gTH_PT(2,nodes_inside_plume);
% % r_nodes_inside_plume     = gTH_PT(3,nodes_inside_plume);
% % theta_plume_center_rot   = repmat(GCOORD_SPH_plume_center_rot(1),1,size(theta_nodes_inside_plume,2));
% % phi_plume_center_rot     = repmat(GCOORD_SPH_plume_center_rot(2),1,size(phi_nodes_inside_plume,2));
% % u_plume = u_max * (1 - ( (r_nodes_inside_plume .* (theta_nodes_inside_plume - theta_plume_center_rot)/r_plume).^2 + ...
% %                          (r_nodes_inside_plume .* (phi_nodes_inside_plume   - phi_plume_center_rot)/r_plume).^2) );
% % 
% % scatter3(MESH.GCOORD_SPH(1,nodes_inside_plume), ...
% %          MESH.GCOORD_SPH(2,nodes_inside_plume), ...
% %          (u_plume/10 + MESH.GCOORD_SPH(3,nodes_inside_plume))/1000,'MarkerEdgeColor','k','MarkerFaceColor',[0.8 0.8 0.8])


% Compute bottom disc in Cartesian coordinates
GCOORD_plume_bottom_center_rot  = spherical2cartesian([GCOORD_SPH_plume_center_rot(1);  GCOORD_SPH_plume_center_rot(2); r_plume_bottom]); 
x_plume_center_rot              = GCOORD_plume_bottom_center_rot(1);
y_plume_center_rot              = GCOORD_plume_bottom_center_rot(2);
z_plume_center_rot              = GCOORD_plume_bottom_center_rot(3);
alpha                           = (0:5:355)*pi/180;
x_plume                         = zeros(1,size(alpha,2));
y_plume                         = zeros(1,size(alpha,2));
for i = 1:size(alpha,2)
   x_plume(i) = x_plume_center_rot + r_plume * cos(alpha(i));
   y_plume(i) = y_plume_center_rot + r_plume * sin(alpha(i));
end
theta_plume_center_rot          = GCOORD_SPH_plume_center_rot(1);
phi_plume_center_rot            = GCOORD_SPH_plume_center_rot(2);
RR = [cos(phi_plume_center_rot)*cos(theta_plume_center_rot)  -sin(phi_plume_center_rot)   cos(phi_plume_center_rot)*sin(theta_plume_center_rot) ; ...
      sin(phi_plume_center_rot)*cos(theta_plume_center_rot)   cos(phi_plume_center_rot)   sin(phi_plume_center_rot)*sin(theta_plume_center_rot) ; ...
                -sin(theta_plume_center_rot)                              0                             cos(theta_plume_center_rot)            ];

GCOORD_plume_bottom             = [x_plume;y_plume;repmat(z_plume_center_rot,1,size(alpha,2))];
GCOORD_plume_bottom_trans       = GCOORD_plume_bottom - repmat(GCOORD_plume_bottom_center_rot,1,size(GCOORD_plume_bottom,2));
GCOORD_plume_bottom_trans_rot   = RR * GCOORD_plume_bottom_trans;
GCOORD_plume_bottom_rot         = GCOORD_plume_bottom_trans_rot + repmat(GCOORD_plume_bottom_center_rot,1,size(GCOORD_plume_bottom,2));

% Compute top disc in Cartesian coordinates
GCOORD_plume_top_center_rot     = spherical2cartesian([GCOORD_SPH_plume_center_rot(1); GCOORD_SPH_plume_center_rot(2); r_plume_top]);
x_plume_center_rot              = GCOORD_plume_top_center_rot(1);
y_plume_center_rot              = GCOORD_plume_top_center_rot(2);
z_plume_center_rot              = GCOORD_plume_top_center_rot(3);
x_plume                         = zeros(1,size(alpha,2));
y_plume                         = zeros(1,size(alpha,2));
for i = 1:size(alpha,2)
   x_plume(i) = x_plume_center_rot + r_plume * cos(alpha(i));
   y_plume(i) = y_plume_center_rot + r_plume * sin(alpha(i));
end
GCOORD_plume_top                = [x_plume;y_plume;repmat(z_plume_center_rot,1,size(alpha,2))];
GCOORD_plume_top_trans          = GCOORD_plume_top - repmat(GCOORD_plume_top_center_rot,1,size(GCOORD_plume_top,2));
GCOORD_plume_top_trans_rot      = RR * GCOORD_plume_top_trans;
GCOORD_plume_top_rot            = GCOORD_plume_top_trans_rot + repmat(GCOORD_plume_top_center_rot,1,size(GCOORD_plume_top,2));

% Search for points inside the plume
GCOORD_plume                    = [GCOORD_plume_bottom_rot GCOORD_plume_top_rot];
EL2NOD_plume_cart               = delaunay(GCOORD_plume');
EL2NOD_plume_cart               = uint32(EL2NOD_plume_cart');
[els2,~,~]                      = tsearch2(GCOORD_plume,EL2NOD_plume_cart,MESH.GCOORD,[],[]);
nodes_inside_plume_cart         = els2 > 0;

% Compute velocities for those points inside the plume
GCOORD_inside_plume             = MESH.GCOORD(:,nodes_inside_plume_cart);
GCOORD_plume_mid_center_rot     = spherical2cartesian([GCOORD_SPH_plume_center_rot(1); GCOORD_SPH_plume_center_rot(2); (r_plume_bottom + r_plume_top)/2]); 
GCOORD_inside_plume_trans       = GCOORD_inside_plume - repmat(GCOORD_plume_mid_center_rot,1,size(GCOORD_inside_plume,2));
GCOORD_inside_plume_trans_unrot = RR' * GCOORD_inside_plume_trans;

Uz_plume_cart = u_max * (1 - ( ((GCOORD_inside_plume_trans_unrot(1,:))/r_plume).^2 + ...
                               ((GCOORD_inside_plume_trans_unrot(2,:))/r_plume).^2) );

% Compute Temperature inside the plume following a gaussian-shaped radial temperature profile
rho                  = sqrt(GCOORD_inside_plume_trans_unrot(1,:).^2 + ...
                           (GCOORD_inside_plume_trans_unrot(2,:)).^2); % radial distance to the axis in the plume center
sigma                = 35;
dT                   = T_plume * ones(sum(nodes_inside_plume_cart),1) - T(nodes_inside_plume_cart);
w_gauss              = exp( -(rho.^2)./(2*sigma^2) ); % gaussian weighting fct to smooth plume temperature
T_nodes_inside_plume = T(nodes_inside_plume_cart) + w_gauss'.*dT;

% Output data
nodes_inside_plume = nodes_inside_plume_cart;
u_plume            = Uz_plume_cart;


figure(348);clf
hold on
axis equal
view(142.5,30)
grid on
% % Plot the shell
% lightGrey           = 0.90*[1 1 1]; % colour for the shell boundaries
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
% scatter3(GCOORD_plume_bottom(1,:), ...
%          GCOORD_plume_bottom(2,:), ...
%          GCOORD_plume_bottom(3,:),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0.8])
% scatter3(GCOORD_plume_bottom_rot(1,:), ...
%          GCOORD_plume_bottom_rot(2,:), ...
%          GCOORD_plume_bottom_rot(3,:),'MarkerEdgeColor','k','MarkerFaceColor',[0.8 0 0])
% scatter3(GCOORD_plume_top(1,:), ...
%          GCOORD_plume_top(2,:), ...
%          GCOORD_plume_top(3,:),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0.8])
% scatter3(GCOORD_plume_top_rot(1,:), ...
%          GCOORD_plume_top_rot(2,:), ...
%          GCOORD_plume_top_rot(3,:),'MarkerEdgeColor','k','MarkerFaceColor',[0.8 0 0])
tetramesh(EL2NOD_plume_cart(1:4,:)',GCOORD_plume','FaceColor',[1 0 0],'FaceAlpha',0.05)
scatter3(MESH.GCOORD(1,nodes_inside_plume_cart), ...
         MESH.GCOORD(2,nodes_inside_plume_cart), ...
         MESH.GCOORD(3,nodes_inside_plume_cart),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1])
     
figure(349);clf
hold on
axis equal
view(142.5,30)
grid on
scatter3(GCOORD_inside_plume_trans_unrot(1,:), ...
         GCOORD_inside_plume_trans_unrot(2,:), ...
         GCOORD_inside_plume_trans_unrot(3,:),'MarkerEdgeColor','k','MarkerFaceColor',[0.8 0 0])
scatter3(GCOORD_inside_plume_trans_unrot(1,:), ...
         GCOORD_inside_plume_trans_unrot(2,:), ...
         GCOORD_inside_plume_trans_unrot(3,:) + Uz_plume_cart,'MarkerEdgeColor','k','MarkerFaceColor',[0.8 0.8 0.8])

figure(350);clf
hold on
axis equal
view(142.5,30)
grid on
scatter3(GCOORD_inside_plume_trans_unrot(1,:), ...
         GCOORD_inside_plume_trans_unrot(2,:), ...
         T_nodes_inside_plume,'MarkerEdgeColor','k','MarkerFaceColor',[0.8 0 0])

end % END OF SUBFUNCTION inside_plume_disc

% #########################################################################

function [velEN,velXYZ] = platevelocity(xyz,platenum,PLATE_VELMODEL)
  % based on veloNE fortran subroutine 12 Mar98 (wjm)
  % The logic is that we construct a vector xyzN pointing towards the North (on
  % the sphere) from the point of interest, and a vector xyzE pointing in
  % the East direction at this point.  (The xyzE vector is constructed at
  % the equator -- it points in the same direction as E at the latitude of
  % the point.) Then the dot product of the Cartesian angular velocity with
  % these direction vectors gives the vector components velN and velE.
  %
  % first find lat,lon of the point (used in logic to find velE and velN componants)
  scale = sum( xyz.^2,1); % dot product -- unrolled loop for l^2 IS THIS LENGTH zero somewhere?
  xyz(:,scale>eps) = xyz(:,scale>eps) ./ repmat( sqrt( scale(scale>eps) ) ,[3 1]); % scale them all by l

  [xyzlatlong] = vangV(xyz);
  % now make a vector xyzN that is 90deg N of the point (wrapping over the pole OK for this logic)
  [xyzN]     = vsetV(xyzlatlong(1)+90.0,xyzlatlong(2));
  % now make a vector pointing east at this point (compute at lat=equator for easier calculation)
  [xyzE]     = vsetV(0.0,xyzlatlong(2)+90.0); 
  % compute angular velocity at the point (cartesian coordinates)
  [velXYZ]   = cross(PLATE_VELMODEL(platenum).Eular_W,xyz);
  [velXYZ]   = velXYZ * 111.111; % conversion from deg/Myr to km/Myr==mm/yr
  velEN(1,:) = dot(xyzE,velXYZ); % velE component
  velEN(2,:) = dot(xyzN,velXYZ); % velN component
end % END OF SUBFUNCTION platevelocity

% #########################################################################

function     [platenum,plateabr] = whichplateV(xyz,PLATES)
%       based on logic in subroutine qplate ( aVlat, along, plate, nplate )
% c              wjm  24jan93   mod. 16june99
% c      CAUTION: No plate can be more than "one hemisphere" across. 
% c        (logic assumes "shortest path" from test-point to boundary)
% 
%     structure for Plate Boundary Table
%     PLATES(iplate).name           = plateFileHeader;
%     PLATES(iplate).plateabr       = plateHeader{iplate}(3:4);
%     PLATES(iplate).nplate         = nplate;
%     PLATES(iplate).platebndlatlon = platebnd{iplate};
%     PLATES(iplate).platebndXYZ    = platebnd{iplate) converted to XYZ;
%
%     xyz input is already in direction-cosine xyz form (instead of lat,lon form)
%     % LOOP OVER ALL PLATE BOUNDARY POINTS TO DO LINE INTEGRAL
%     % ONLY EXIT LOOP IF A POINT IS ON THE BOUNDARY POINT (then in this plate)
%     % 
%     % c--------------------------------------------------
%     % c  A  and  B  are unit vector 90 degrees from  xyz point
%     % c  (found with unit-cross of two boundary points with xyz)
%     % c  (We might prefer to have unit vectors 90 deg from 'xyz' exactly
%     % c  IN THE DIRECTION of the two boundary points from xyz, but being 
%     % c  90 degree shifted by the cross product makes no difference.)
%     % c
%     % c  The magnitude of A cross B  (found by dotting with xyz, to which 
%     % c  it is parallel/antiparallel)  gives 'sine' of the angle we want.
%     % c
%     % c  But if delta-ang is more than 90?  Arcsin will be wrong.(quadrant)
%     % c  (A dot B)  will give the 'cosine' of this angle, so with atan2, OK !!
%     % 
%     % c--------------------------------------------------
% 
%===================================
% BEGINNING OF MAIN LOOP OVER PLATES
%===================================
% make sure xyz is a UNIT vector (=direction cosine)
  scale  = sum( xyz.^2,1); % dot product -- unrolled loop for l^2 IS THIS LENGTH zero somewhere?
  xyz(:,scale>eps) = xyz(:,scale>eps) ./ repmat( sqrt( scale(scale>eps) ) ,[3 1]); % scale them all by l

% loop through all the plates
nplate = length(PLATES);
for iplate = 1:nplate
    % c            if there is an 'exit' (point on node or on connecting line)
    % c            "this-plate" will be the answer  (the one tested for first,
    % c            so the order of the plates in 'PLATES.dat' can matter).
    plateabr = PLATES(iplate).plateabr;
    platenum = iplate;
    nBndPnt  = size(PLATES(iplate).platebndXYZ,2);
    
    xyzPnt = repmat(xyz,1,nBndPnt);
    A      = cross(xyzPnt,PLATES(iplate).platebndXYZ);
    scale  = dot(A,A);
    ionBnd = find( scale < 1e-10 );
    if ~isempty(ionBnd)
        if dot(xyz,A(:,ionBnd))>0.99
            break
        end
    end
    scale  = sqrt(scale);
    A      = A./repmat(scale,3,1); % makes all A's unit vectors
    B      = A(:,2:end);
    C      = cross(A(:,1:end-1),B);
    
    sina   = dot(C,xyzPnt(:,1:end-1));
    cosa   = dot(A(:,1:end-1),B);
    delang = atan2(sina,cosa);
    
    enclosed_angle = sum(delang);
    if enclosed_angle >= 6.275 % this is the 2 Pi case, exit with plate number 
        break
    end
% normal exit, return with a plate name
%======================
end; % END OF MAIN LOOP
%======================

end % END OF SUBFUNCTION whichplate

% #########################################################################

%---------------------------------------------
% function [V] = cross(A,B) % replaced by built-in MATLAB cross-product
% %       find cross product  V = A x B
%   V = zeros(size(A));
%   V(1,:) = A(2,:) .* B(3,:)  -  A(3,:) .* B(2,:);
%   V(2,:) = A(3,:) .* B(1,:)  -  A(1,:) .* B(3,:);
%   V(3,:) = A(1,:) .* B(2,:)  -  A(2,:) .* B(1,:);
% end
%---------------------------------------------

function [U,AeqB] = ucross(A,B)
%   find unit  vector perpendicular to A and B, U= unit(A cross B)
%   also return logical vector that is 1 when A=B to machine 10*eps (then cross=0)
  U = zeros(size(A));
  U = cross(A,B); % take cross products -- U vectors pointing perp to A and B
  % now find length of all U to rescale to unit vector
  scale    = sum( U.^2,1); % dot product -- unrolled loop for l^2
  AeqB = scale < 10*eps;  % logical find where scale is zero (i.e. when A=B)
  U(:,~AeqB) = U(:,~AeqB) ./ repmat( sqrt(scale(~AeqB)) ,[3 1]); % scale them all by l
                                % repmat 3 copies of l, 1 for each component
end % END OF SUBFUNCTION ucross

% #########################################################################

function [alatlong] = vangV(a)
%c        convert from direction cosines into  lat,long   26apr99
%c        (vector 'a' does not have to be of unit length.)
  rad2deg =  1./0.0174532925199;
  aa = zeros(1,size(a,2));   % initialize variables
  alat= zeros(1,size(a,2));
  along = zeros(1,size(a,2));
  
  aa(:) = sqrt( a(1,:).*a(1,:) + a(2,:).*a(2,:) );
  alat(:) = 90.0 .* sign( a(3,:) );
  along(:) = 0.0;
  
  if aa == 0
      alat(:)  = 90.0 .* sign( a(3,:) );
      along(:) = 0;
  else
      alat(aa > 0.0)  = atan ( a(3,aa > 0.0) ./ aa(aa > 0.0)  ) .* rad2deg;
      along(aa > 0.0) = atan2( a(2,aa > 0.0), a(1,aa > 0.0) ) .* rad2deg;
  end

%c...(uncomment next lines to go from 0-360 degrees instead of -180 to +180)
  along(along < 0.0 ) = along(along < 0.0 ) + 360.;
  along(along > 359.999) = 0.;
  alatlong = [alat;along];
end % END OF SUBFUNCTION vangV

% #########################################################################

function [a] = vsetV(alat,along)
%c convert  lat,long in degrees into direction cosines  26apr99
  a = zeros(3,max(size(alat)));
  deg2rad =  0.0174532925199;
  a(1,:) = cos(deg2rad*alat(:)) .* cos(deg2rad*along(:));
  a(2,:) = cos(deg2rad*alat(:)) .* sin(deg2rad*along(:));
  a(3,:) = sin(deg2rad*alat(:));
end % END OF SUBFUNCTION vsetV
