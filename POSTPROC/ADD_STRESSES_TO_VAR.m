function ADD_STRESSES_TO_VAR(OPTS)
% Usage: ADD_STRESSES_TO_VAR(OPTS)
%
% Purpose: 
%   Add stresses to output data where stresses were not calculated
%
% Input:
%   OPTS : [structure] : Structure containing info for output file
%
% Output:
%
% JMT Apr 2017


pdir = pwd;
cd('..');
addpath([pwd '/mfiles_M3TET']);
addpath([pwd '/mfiles_SPH']);
cd(pdir);

%==========
% DEFAULTS
%==========
% (1) specify output directory (folder where data is located)
outdir = '/Users/jorge/Tests/SPH/Trash_01';
% (2) 
nelblk = 1400;
% (2) 
dt     = 0;

% Overwrite defaults if different values are specified in structure OPTS
if nargin>0
    if isfield(OPTS,'outdir')
        outdir = OPTS.outdir;
    end
    if isfield(OPTS,'nelblk')
        outdir = OPTS.nelblk;
    end
    if isfield(OPTS,'dt')
        outdir = OPTS.dt;
    end
end

%==========================================================================
% LOAD MESH DATA
%==========================================================================
outdir  = [outdir filesep];
fprintf(' Loading mesh data and settings...');
try
    data  = load([outdir 'Lab01x01_Mesh.mat']);
catch
    error(' Could no open file "%s"\n',[outdir 'MESH.mat']);
end
if isfield(data,'Mesh');
    MESH = data.Mesh;
else
    MESH = data.MESH;
end
clear data
if isfield(MESH,'Phase_id')
    MESH.PhaseID = MESH.Phase_id;
elseif ~isfield(MESH,'PhaseID')
    MESH.PhaseID = 1;
end

%==========================================================================
% LOAD SETTINGS
%==========================================================================
try
    data2  = load([outdir 'Lab01x01_SETTINGS.mat']);
catch
    error(' Could no open file "%s"\n',[outdir 'SETTINGS.mat']);
end
SETTINGS = data2.SETTINGS;
PHYSICS  = data2.PHYSICS;
NUMSCALE = data2.NUMSCALE;
clear data2
SETTINGS.is_elastic = 'no';
fprintf(' done\n');

%==========================================================================
% CALCULATE STRESSES AND ADD THEM TO VAR
%==========================================================================
cd(outdir);
Data_info   = dir('Lab01x01_Data_*');
num_files   = numel(Data_info); % number of Data_info files
g           = abs(PHYSICS.g); % [m/s^2]
Dens_mantle = PHYSICS.Dens;
Dens_water  = 1030; % [kg/m^3]

for i = 1:num_files
    % Store subdomain VAR and time
    load(Data_info(i).name)
    
%     VAR.sigma_n = [];
%     
%     VAR_temp = calc_stresses_3d(VAR,MESH,SETTINGS,PHYSICS,NUMSCALE,dt,nelblk);
%     
%     %======================================================================
%     % MAPPING FROM INTEGRATION POINTS TO NODES - IP CANNOT BE PLOTTED
%     %======================================================================
%     
% % %     % USING SCATTEREDINTERPOLANT BUILT-IN MATLAB FUNCTION
% % %     x_ip        = VAR_temp.x_ip(:);
% % %     y_ip        = VAR_temp.y_ip(:);
% % %     z_ip        = VAR_temp.z_ip(:);
% % %     sigma_n     = VAR_temp.sigma_n(:);
% % %     % ------------------- Interpolation functions
% % %     F_sigma_n = scatteredInterpolant(x_ip,y_ip,z_ip,sigma_n);
% % %     
% % %     % ------------------- Evaluate variables at query points
% % %     sigma_n_q   = F_sigma_n(MESH.GCOORD(1,1:MESH.nnod), ...
% % %                             MESH.GCOORD(2,1:MESH.nnod), ...
% % %                             MESH.GCOORD(3,1:MESH.nnod));
% % %     VAR.sigma_n = sigma_n_q';                  
%     
%     % USING JOERG'S FUNCTION: - much faster
%     %                         - results a bit different than using scatteredInterpolant
%     VAR.sigma_n = ipval_to_nodval_3d(MESH.GCOORD,MESH.EL2NOD{1},VAR_temp.sigma_n);
%     
%     %======================================================================
%     % COMPUTING DYNAMIC TOPOGRAPHY
%     %======================================================================
%     VAR.h_dyn = VAR.sigma_n/((Dens_mantle - Dens_water)*g); % [m]
%     VAR.h_dyn(MESH.PointID ~= 306) = 0; % save dynamic topography only for surface points
    
    %======================================================================
    % COMPUTING MANTLE CONTRIBUTION TO TOPOGRAPHY (ISOSTATIC TOPOGRAPHY)
    %======================================================================
    VAR.h_iso = isostatic_integration(MESH,VAR,PHYSICS);
    
    suffix = ['_' num2str_d(i-1,4)];
    save([outdir 'Lab01x01_Data' suffix],'dt','istep','time','VAR');
end
cd(pdir)

end % END OF FUNCTION ADD_STRESSES_TO_VAR

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function h_iso = isostatic_integration(MESH,VAR,PHYSICS)
% Usage: h_iso = isostatic_integration(MESH,VAR,PHYSICS)
%
% Purpose: 
%   Compute the mantle contribution to topography. First, the thermal
%   subsidence is computed making an isostatic integration. For this, we
%   create a layer with surface points and then we add deeper layers
%   projecting those surface points into the radial direction (so that we
%   have sample columns at each surface point). 
%
% Input:
%   MESH    : [structure] : structure containing FE mesh data
%   VAR     : [structure] : structure containing all major variables
%   PHYSICS : [structure] : structure containing physical properties
%
% Output:
%   h_iso   : [vector]    : mantle contribution to topography
%
% JMT Apr 2017

%==========================================================================
% DEFINE VARIABLES
%==========================================================================

Dens_mantle      = PHYSICS.Dens; % [kg/m^3]
Dens_water       = 1030;         % [kg/m^3]
h_iso            = zeros(MESH.nnod,1);

GCOORD_SPH_surf  = MESH.GCOORD_SPH(:,MESH.PointID == 306); % surface points
nnod_surf        = size(GCOORD_SPH_surf,2);
layer_depth      = 300; % Compensated depth [km] 
step_depth       = 10;  % [km]
layer_step_depth = [zeros(2,nnod_surf); step_depth*ones(1,nnod_surf)];

GCOORD_SPH_this_depth = GCOORD_SPH_surf;
layers                = round(layer_depth/step_depth) + 1;
GCOORD_SPH_sample     = zeros(3,nnod_surf,layers);
GCOORD_sample         = zeros(3,nnod_surf,layers);
Dens_q                = zeros(1,nnod_surf,layers);

% USING SCATTEREDINTERPOLANT BUILT-IN MATLAB FUNCTION
x        = MESH.GCOORD(1,:)';
y        = MESH.GCOORD(2,:)';
z        = MESH.GCOORD(3,:)';
Dens     = VAR.Dens(:);
% ------------------- Interpolation functions
F_Dens = scatteredInterpolant(x,y,z,Dens);

for i = 1:layers
    GCOORD_SPH_sample(:,:,i) = GCOORD_SPH_this_depth;
    GCOORD_sample(:,:,i)     = spherical2cartesian(GCOORD_SPH_sample(:,:,i));
    % ------------------- Evaluate variables at query points
    Dens_q(1,:,i)            = F_Dens(GCOORD_sample(1,:,i), ...
                                      GCOORD_sample(2,:,i), ...
                                      GCOORD_sample(3,:,i));
    GCOORD_SPH_this_depth    = GCOORD_SPH_this_depth - layer_step_depth;
end

Dens_temp  = squeeze(Dens_q)';

%==========================================================================
% RIEMANN INTEGRAL TO COMPUTE THERMAL SUBSIDENCE (pick one)
%==========================================================================
% h_iso_surf_left  = ...
%     sum((Dens_temp(1:layers-1,:) - Dens_mantle)*step_depth)/(Dens_mantle - Dens_water); % [km]
h_iso_surf_right = ...
    sum((Dens_temp(2:layers,:) - Dens_mantle)*step_depth)/(Dens_mantle - Dens_water); % [km]

%==========================================================================
% MANTLE CONTRIBUTION TO TOPOGRAPHY
%==========================================================================
% Change the sign because positive subsidence = negative topography
h_iso(MESH.PointID == 306) = -h_iso_surf_right; 

end % END OF SUBFUNCTION isostatic_integration