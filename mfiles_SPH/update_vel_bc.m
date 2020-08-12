function [UBC,VAR] = update_vel_bc(MESH,COMM,PHYSICS,VAR,UBC,model_time)
% Usage: [UBC,VAR] = update_vel_bc(MESH,COMM,PHYSICS,VAR,UBC,model_time)
%
% Purpose: Update velocity boundary conditions from GPlates files. The
%          velocities are linearly interpolated and weighted depending of
%          how close is model_time to the ends of the time interval given
%          by the GPlates files (usually 1 Myr) 
%
% Input:
%   MESH       : [structure] : structure containing FE mesh data
%   COMM       : [structure] : inter-subdomain communication data
%   PHYSICS    : [structure] : structure containing physical properties
%   VAR        : [structure] : structure containing all major variables
%   UBC        : [structure] : structure containing velocity boundary conditions
%   model_time : [scalar]    : time in the running model
%
% Output:
%   UBC        : [structure] : structure containing velocity boundary conditions
%   VAR        : [structure] : structure containing all major variables
%
% JMT Feb 2017

fprintf('\n Updating GPlates velocity BCs...');
%==========================================================================
% LOAD VELOCITY BCs AND TIMES
%==========================================================================
V_BCs           = PHYSICS.V_BCs;
time_VBCs_files = PHYSICS.time_VBCs_files;
model_time      = max(time_VBCs_files) - model_time; % write model time in terms of the time for the plate reconstruction

%==========================================================================
% INTERPOLATE VELOCITIES
%==========================================================================
i_files = find(abs(time_VBCs_files- model_time) < 1); % check in which interval of "time_files" is "model_time"
if size(i_files,1) == 1
    if ~(i_files == 1 || i_files == size(V_BCs,2))
        weight  = (model_time - time_VBCs_files(i_files)) / (time_VBCs_files(i_files + 1) - time_VBCs_files(i_files));
        Uth     = V_BCs(i_files).Uth_Uph(1,:) * (1 - weight) + V_BCs(i_files + 1).Uth_Uph(1,:) * weight;
        Uph     = V_BCs(i_files).Uth_Uph(2,:) * (1 - weight) + V_BCs(i_files + 1).Uth_Uph(2,:) * weight;
        Uth_Uph = [Uth; Uph];
    else
        Uth_Uph = V_BCs(i_files(1)).Uth_Uph;
    end
else
    weight  = (model_time - time_VBCs_files(i_files(1))) / (time_VBCs_files(i_files(2)) - time_VBCs_files(i_files(1)));
    Uth     = V_BCs(i_files(1)).Uth_Uph(1,:) * (1 - weight) + V_BCs(i_files(2)).Uth_Uph(1,:) * weight;
    Uph     = V_BCs(i_files(1)).Uth_Uph(2,:) * (1 - weight) + V_BCs(i_files(2)).Uth_Uph(2,:) * weight;
    Uth_Uph = [Uth; Uph];
end

%==========================================================================
% UPDATE THE VELOCITY BCs WITH INTERPOLATED VELOCITIES
%==========================================================================
% Make sure that all domain boundary indices, on which boundary conditions will be imposed, have been located in the mesh
indx    = [301 306];
DBnods  = prepare_DBnods(MESH,COMM,indx);
DBindx  = MESH.PointID(DBnods);
nDBnods = length(DBnods);

% Allocate arrays
iUbc      = 0;                  % counter for velocity boundary conditions
iUfix     = zeros(1,3*nDBnods); % list of prescribed velocity dofs
vUfix     = zeros(1,3*nDBnods); % list of the prescribed velocity values
iUth      = 0;                  % counter for reading Uth_Uph matrix
iUph      = 0;                  % counter for reading Uth_Uph matrix

% Rotate GPlates velocities (actual Earth) into the spherical mesh frame (where the refined region is centered at theta = 90, phi = 90)
sph             = PHYSICS.V_BCs(1).latlon;
sph(1,:)        = 90 - sph(1,:); % transform latitude to colatitude
sph             = (pi/180)*sph;  % transform degrees to radians
sph             = [sph; MESH.r_surf*ones(1,size(sph,2))];
cart            = spherical2cartesian(sph);
TT              = make_rotation_matrix(cart); % transformation matrix between spherical and Cartesian coordinates
u_sph           = [Uth_Uph; zeros(1,size(Uth_Uph,2))];
u_cart          = TT * u_sph(:);
u_cart_rot      = (MESH.RR2 * MESH.RR1)' * reshape(u_cart,3,[]); % rotate velocities
cart_rot        = (MESH.RR2 * MESH.RR1)' * cart;  % rotate coordinates
TT_rot          = make_rotation_matrix(cart_rot); % transformation matrix between spherical and Cartesian rotated coordinates
u_sph_rot       = TT_rot'* u_cart_rot(:);
u_sph_rot       = reshape(u_sph_rot,3,[]);
u_sph_rot(3,:)  = [];
Uth_Uph         = u_sph_rot;

% Loop over all nodes on domain boundary
for ii=1:nDBnods
    node   = DBnods(ii); % global node number (= global temperature dof)
    ndf_th = 3*node-2;   % global 1st velocity dof (colatitude (theta) direction)
    ndf_ph = 3*node-1;   % global 2nd velocity dof (longitude (phi) direction)
    ndf_r  = 3*node;     % global 3rd velocity dof (radial (r) direction)
    
    switch DBindx(ii)
        case 301 % CORE-MANTLE-BOUNDARY
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
            % Velocity BC
            iUth        = iUth+1;
            iUbc        = iUbc+1;
            iUfix(iUbc) = ndf_th; % dof colatitude-comp of Vel (number)
            vUfix(iUbc) = Uth_Uph(1,iUth); % colatitude-comp of Vel (value)
            
            iUph        = iUph+1;
            iUbc        = iUbc+1;
            iUfix(iUbc) = ndf_ph; % dof longitude-comp of Vel (number)
            vUfix(iUbc) = Uth_Uph(2,iUph); % longitude-comp of Vel (value)
            
            iUbc        = iUbc+1;
            iUfix(iUbc) = ndf_r; % this is the rotated
            vUfix(iUbc) = 0;     % tangent motion to the surface
    end
end

% shrink arrays to correct size
iUfix(iUbc+1:end)       = [];
vUfix(iUbc+1:end)       = [];

if ~isempty(UBC.PLUME)
    % LOAD VELOCITY BCS FOR A PLUME COMPUTED IN bc_sph
    % ================================================
    iUbc_plume  = UBC.PLUME.iUbc_plume;
    iUfix_plume = UBC.PLUME.iUfix_plume;
    vUfix_plume = UBC.PLUME.vUfix_plume;
else
    iUbc_plume  = 0;
    iUfix_plume = [];
    vUfix_plume = [];
    UBC.PLUME   = [];
end

% Update boundary condition structures
UBC = struct('nfix'      , iUbc + iUbc_plume,...
             'ifix'      , [iUfix iUfix_plume],...
             'vfix'      , [vUfix vUfix_plume],...
             'PlateInfo' , PHYSICS.V_BCs(i_files(1)).PlateID,... % updated PlateID
             'PLUME'     , UBC.PLUME);
VAR.Plate                      = zeros(MESH.nnod,1);
VAR.Plate(MESH.PointID == 306) = UBC.PlateInfo;

fprintf(' done.\n\n');

end % END OF FUNCTION update_vel_bc

% #########################################################################
%                              SUB-FUNCTIONS
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