function [PHYSICS,SETTINGS,NUMSCALE] = model_parameters_sph_south_atlantic_plume(plume_flux,colat_lon_plume_center)
% Usage: [PHYSICS,SETTINGS,NUMSCALE] = model_parameters_sph_south_atlantic_plume(plume_flux,colat_lon_plume_center)
%
% Purpose: Define model parameters and physical properties
%
% Input:
%   plume_flux             : [scalar] : plume flux (km^3/yr)
%   colat_lon_plume_center : [vector] : % center of the plume
%                                       (colatitude,longitude) [degrees]
%
% Output:
%   SETTINGS : [structure] : model parameters
%   PHYSICS  : [structure] : physical properties
%   NUMSCALE : [structure] : numerical scaling parameters
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Mar 2012
% JH Feb 2016 : updated to newest version
% JMT Nov 2018 : updated to select Meshes and GPlates velocity model 
%

%==========================================================================
% CODE SETUP AND DEFINITION OF MATERIAL PROPERTIES
%==========================================================================

% NUMSCALE defines the units in which the code operates. If not indicated
% otherwise all parameters must by provided in these units.
% Furthermore, NUMSCALE.Dens0 defines the reference density, meaning that
% PHYSICS.g*NUMSCALE.Dens0 is subtracted as lithostatic pressure.
NUMSCALE = numerical_eq_scaling('mantle'); % "mantle"  or "SI


%==========================================================================
% LOAD DEFAULTS FOR SOLVERS AND ADVECTION SCHEME
SETTINGS = default_numerics_M3TET;
%==========================================================================


%==========================================================================
% START NEW CALCULATION OR RESTART?
SETTINGS.restart       = 0;
  % 0 --> start new calculation
  % 1 --> continue a calculation
SETTINGS.restart_file  = 'last';
  % Which file shall be used for a restart? Allowed values are:
  % 'last'   --> uses the most recent restart file
  % a number --> the number of the Data_xxxx.mat file to be loaded
  %              (e.g. restart_file = 17 will load "Data_0017.mat")
SETTINGS.restart_setup = 'from_file';
  % Which settings shall be used for a restart? Allowed values are:
  % 'from_file'   --> PHYSICS and SETTINGS will be loaded from file
  %                   (everything defined in this function will be
  %                   overwritten)
  % 'new_defined' --> Only variables will be loaded from file Data_xxxx.mat
  %                   PHYSICS and SETTINGS are defined below and NOT loaded
  %                   from file. This allows you to change some settings.
  %                   NOTE: Changing xmin/xmax/zmin/zmax will have no effect
  %                         since the old mesh is loaded!!!!
%==========================================================================


%==========================================================================
% OUTPUT FOLDER, OUTPUT INTERVAL AND SIMULATION TIME
path2data            = data_storage_3d();
  % see m-file "data_storage_3d.m" to add location on new computer
plume = strcat(num2str(plume_flux),'_',num2str(colat_lon_plume_center(1)),'_',num2str(colat_lon_plume_center(2)));
if plume_flux == 0
    SETTINGS.outdir = [path2data 'SPH/NoPlume'];
else
    switch plume
        case '10_118_354'  % L1F10
            SETTINGS.outdir = [path2data 'SPH/Plume_flux_10_Plume_position_118_354'];
        case '15_118_354'  % L1F15
            SETTINGS.outdir = [path2data 'SPH/Plume_flux_15_Plume_position_118_354'];
        case '20_118_354'  % L1F20
            SETTINGS.outdir = [path2data 'SPH/Plume_flux_20_Plume_position_118_354'];
        case '5_118_354'   % L1F05
            SETTINGS.outdir = [path2data 'SPH/Plume_flux_05_Plume_position_118_354'];
        case '10_121_351'  % L2F10
            SETTINGS.outdir = [path2data 'SPH/Plume_flux_10_Plume_position_121_351'];
        case '7.5_127_348' % L4F17.5
            SETTINGS.outdir = [path2data 'SPH/Plume_flux_7_5_Plume_position_127_348'];
        case '15_124_351'  % L3F15
            SETTINGS.outdir = [path2data 'SPH/Plume_flux_15_Plume_position_124_351'];
        case '15_121_351'  % L2F15
            SETTINGS.outdir = [path2data 'SPH/Plume_flux_15_Plume_position_121_351'];
        case '20_121_351'  % L2F20
            SETTINGS.outdir = [path2data 'SPH/Plume_flux_20_Plume_position_121_351'];
        case '15_120_0'    % L1bF15 (plume beneath Africa)
            SETTINGS.outdir = [path2data 'SPH/Plume_flux_15_Plume_position_120_0'];
        otherwise
            error('create new case')
    end
end
  % full path to output folder
SETTINGS.save_period = 0.1;
  % [in units of NUMSCALE.t0]; save/plot interval
% Simulation stops after either endtime or maximum number of time steps
% has been reached.
SETTINGS.endtime     = 30.1; %30.1;
  % [in units of NUMSCALE.t0]; duration of simulation
SETTINGS.nstep       = 1e8;
  % max number of time steps
%==========================================================================


%==========================================================================
% FIGURE PLOTTING
SETTINGS.save_figs   = 1;     % 1--> saves figures in output directory
SETTINGS.fig_type    = 'png'; % 'png' or 'fig'; format of saved figure
SETTINGS.show_figs   = 1;     % 0--> figures not shown on screen (faster)
                              % 1--> figures shown on screen
%==========================================================================


%==========================================================================
% TIME STEP CONSTRAINTS
SETTINGS.dtmax   = 0.1; % 0.05
  % [in units of NUMSCALE.t0]; maximum (user-defined) time step
  % set dtmax==0 to disable
SETTINGS.dt_adv_limit = 0.005;
  % [in units of NUMSCALE.t0]; max advection distance during one time step
  % It limits the time step dt so that nothing advects more than this
  % fraction of Lx, Ly and Lz
  % Set to zero to limit advection by the size of an element (then 
  % nothing advects further than an element length; can be slow...)
%==========================================================================


%==========================================================================
% MESH FILE, ELEMENT TYPE, EDGE TYPE AND DOMAIN SIZE
SETTINGS.meshfile = 'SETUP_TEST/MESHES/EarthShell_n37k';
  % Regular meshes:
  %     EarthShell_reg_n2k  -> L0 = 2000 km 
  %     EarthShell_reg_n4k  -> L0 = 1500 km 
  %     EarthShell_reg_n5k  -> L0 = 1400 km 
  %     EarthShell_reg_n7k  -> L0 = 1300 km 
  %     EarthShell_reg_n9k  -> L0 = 1200 km 
  %     EarthShell_reg_n11k -> L0 = 1100 km 
  %     EarthShell_reg_n13k -> L0 = 1000 km 
  %     EarthShell_reg_n17k -> L0 =  900 km 
  %     EarthShell_reg_n23k -> L0 =  800 km 
  %     EarthShell_reg_n35k -> L0 =  700 km 
  %     EarthShell_reg_n61k -> L0 =  600 km 
  %     EarthShell_reg_n93k -> L0 =  500 km 
  %
  % Meshes with embedded region:
  %                     ------------------------------------------------------------------------------------------
  %                     | Region     | % nodes | % elements | element lenght | width(N-S) | width(E-W) |  depth  |
  %   ------------------------------------------------------------------------------------------------------------
  %   |                 | Coarse     |   7.4   |     6.3    |    2000 km     |      -     |      -     |    -    |
  %   | EarthShell_n26k | Transition |  36.2   |    37.6    | 2000 to 150 km |   8000 km  |   8000 km  | 2900 km |
  %   |    (squared)    | Refined    |  56.4   |    56.2    |     150 km     |   3333 km  |   3333 km  |  300 km |
  %   ------------------------------------------------------------------------------------------------------------
  %   |                 | Coarse     |   5.2   |     4.4    |    2000 km     |      -     |      -     |    -    |
  %   | EarthShell_n37k | Transition |  33.4   |    29.7    | 2000 to 100 km |   8000 km  |   5600 km  | 2900 km |
  %   |  (rectangular)  | Refined    |  64.8   |    65.9    |     100 km     |   4200 km  |   1800 km  |  300 km |
  %   ------------------------------------------------------------------------------------------------------------
  %   |                 | Coarse     |   2.9   |     2.4    |    2000 km     |      -     |      -     |    -    |
  %   | EarthShell_n64k | Transition |  27.6   |    28.2    | 2000 to 100 km |   9600 km  |   6800 km  | 2900 km |
  %   |  (rectangular)  | Refined    |  69.5   |    69.4    |     100 km     |   5000 km  |   2200 km  |  300 km |
  %   ------------------------------------------------------------------------------------------------------------
  %   |                 | Coarse     |   2.9   |     2.4    |    2000 km     |      -     |      -     |    -    |
  %   | EarthShell_n66k | Transition |  28.7   |    29.5    | 2000 to 100 km |   8000 km  |   8000 km  | 2900 km |
  %   |    (squared)    | Refined    |  68.3   |    68.1    |     100 km     |   3333 km  |   3333 km  |  300 km |
  %   ------------------------------------------------------------------------------------------------------------
  %
  % Specify path to and the name of file containing the mesh data in format
  % MESH.GCOORD, MESH.EL2NOD, MESH.PointID, (and optionally MESH.PhaseID)

SETTINGS.nmg = 3;
  % Define number of multigrid levels (nmg<2 will disable multigrid)
  % CAUTION: With each level the number of nodes increases by a factor of
  %          about 7-8: If the loaded mesh has 100k nodes, the final mesh
  %          with nmg=2 will have ~750k nodes and 5.6 million with nmg=3

SETTINGS.GPlates_Vel_BCs_file = 'Vel_BCs_n37k_nmg3'; % Name of the file with GPlates velocity data

SETTINGS.element_type = 'quadratic'; % 'cubic' need to be coded!!
  % Define which kind of element will be used
  % 'quadratic' --> 10 nodes/element
  % 'cubic'     --> 20 nodes/element

SETTINGS.edge_type = 'curved';
  % Define which kind of edges will be used for the elements of the mesh
  % 'curved'   --> curved edges
  % 'straight' --> straight edges (not useful for spherical meshes)

SETTINGS.jacobian = 'double';
  % Define which method will be used for assembling the stiffness matrices
  % 'standard' --> standard method 
  % 'double'   --> double jacobian approach (faster and more accurate than
  %                standard method)
  
SETTINGS.theta_cone = 45; % half aperture of the cone. theta_cone is 
                          % measured from the positive Z axis (North Pole)
                          % towards the South Pole (0? < theta_cone < 90?)                          
  
PHYSICS.xmin = [];  % start of x-axis
PHYSICS.xmax = [];  % end of x-axis
PHYSICS.ymin = [];  % start of y-axis
PHYSICS.ymax = [];  % end of y-axis
PHYSICS.zmin = [];  % bottom of domain
PHYSICS.zmax = [];  % top of domain
  % [in units of NUMSCALE.L0]
  % (physical size of domain will be PHYSICS.xmax*NUMSCALE.L0, etc)
%==========================================================================


%==========================================================================
% INFO ON TECTONIC PLATES, VELOCITY MODEL AND PLUME
SETTINGS.plate_model = 'GPlates';
  % WJM     --> Plate definition (William Jason Morgan) and 
  %             Hotspot velocity model (Morgan & Phipps Morgan, 2007)
  % GPlates --> Plate kinematic reconstructions (Gurnis et al., 2012)
  % Debug   --> Simple case for debugging
if strcmp(SETTINGS.plate_model,'GPlates')
  if plume_flux == 0
      SETTINGS.plume = 'no';  % Only coded for SETTINGS.plate_model = 'GPlates'
  else
      SETTINGS.plume = 'yes'; % Only coded for SETTINGS.plate_model = 'GPlates'
  end
end
SETTINGS.cratons_contours = '8M';
  % '3'  --> 3 countours to define the shape of each craton. (OLD VERSION)
  % '8'  --> 8 countours to define the shape of each craton and diffusion 
  %          is run for 10 Myr before start the simulation to smooth even  
  %          more the shape of the cratons. (OLD VERSION)
  % '8M' --> 8 countours to define the shape of each craton. Congo and Sao
  %          Francisco cratons are merged in one contour. Kalahari and Rio 
  %          de la Plata cratons are merged in one contour. Diffusion is
  %          run for 10 Myr before start the simulation to smooth even more
  %          the shape of the cratons
SETTINGS.del_lith_thickness = 0; % difference of lithospheric thickness (in km)
  % from latitude 118° to 145° (~ 3000 km). Allowed values: -60 km to 60 km
  % del_lith_thickness = 0 --> No slope
  % del_lith_thickness > 0 --> Makes the lithopshere thicker towards the South 
  % del_lith_thickness < 0 --> Makes the lithopshere thinner towards the South
if SETTINGS.del_lith_thickness < -60 || SETTINGS.del_lith_thickness > 60
   error ('del_lith_thickness must be between -60 and 60') 
end
switch SETTINGS.plate_model
    case 'WJM'
        SETTINGS.useNataf = 0;
          % 0 --> plate motion prescribed everywhere
          % 1 --> only Nataf craton motion prescribed

        % READ PLATE BOUNDARY INFO
        addpath('./SETUP_TEST/PLATES');
        PlateFileName          = './SETUP_TEST/PLATES/PLATES.dat';
        PHYSICS.PLATES         = readPlateBndInfo(PlateFileName);
          % structure for Plate Boundary Table
          % PLATES(1:nplate) = struct('name'          ,' ',...
          %                           'plateabr'      ,' ',...
          %                           'nplate'        ,nplate,...
          %                           'platebndlatlon',zeros(2),...
          %                           'platebndXYZ',zeros(3));

        % READ PLATE VELOCITY MODEL
        VelModelFileName       = './SETUP_TEST/PLATES/VELOCITY_MODEL.txt';
        PHYSICS.PLATE_VELMODEL = readPlateVelModel(VelModelFileName);
          % structure for plate velocity Eular Pole table
          % PLATE_VELMODEL(1:nplate) = struct('name'      ,' ',...
          %                           'platenum',nplate,...
          %                           'plateabr'      ,' ',...
          %                           'nplate'        ,nplate,...
          %                           'EularLatLonAngrate',zeros(3),...
          %                           'Eular_W',zeros(3),...
          %                           'Eqrate',zeros(1));

        NatafTectRegFileName   = './SETUP_TEST/PLATES/NATAF2BY2.txt';
        PHYSICS.NATAF_CRATON   = readNatafCratonInfo(NatafTectRegFileName);

    case 'GPlates'
        % READ GPLATES VELOCITY MODEL
        % ===========================
        pdir = pwd;
        cd SETUP_TEST/GPLATES/VEL_BCs;
        load(SETTINGS.GPlates_Vel_BCs_file)
        PHYSICS.V_BCs           = V_BCs;
        PHYSICS.time_VBCs_files = time_VBCs_files;

        % READ CRATONS INFO
        % =================
        switch SETTINGS.cratons_contours
            case '3'
                cd ..; cd CRATONS_3_contours
            case '8'
                cd ..; cd CRATONS_8_contours;
            case '8M'
                cd ..; cd CRATONS_MERGED_8_contours;
            otherwise
                error('SETTINGS.cratons_contours must be ''3'', ''8'' or ''8M''')
        end
        %------------------------------------------------------------------
        % Check Cratons file format and number of files
        %------------------------------------------------------------------
        Cratons_info = dir();                % get info about the files in this directory
        sample_file  = Cratons_info(4).name; % select 1st craton file to check the extension
        if strcmp(sample_file(end-2:end),'.xy')
            system('./xy2txt.sh'); % change .xy extension to .txt extension
        elseif strcmp(sample_file(end-3:end),'.txt')
            % Cratons files are already in .txt format
        else
            error('wrong extension for CRATONS files. It should be .xy or .txt')
        end
        Cratons_info = dir('*.txt');
        num_files    = numel(Cratons_info); % number of Cratons files
        
        %------------------------------------------------------------------
        % Create structure for CRATONS and load the files
        %------------------------------------------------------------------
        CRATONS(1:num_files) = struct('latlon',0);
        fprintf('\n Loading Cratons info...');
        for i = 1:num_files
            fid               = fopen(Cratons_info(i).name,'r'); % open input file
            tmp               = fscanf(fid, '%f');
            CRATONS(i).latlon = reshape(tmp,2,[])';
            fclose(fid); % close input file
        end
        fprintf(' done.\n\n');
        cd(pdir);
        PHYSICS.CRATONS = CRATONS;
        clear Cratons_info i numfiles pdir sample_file
        
        switch SETTINGS.plume
            case 'yes'
                % SET INFO FOR TRISTAN DA CUNHA PLUME (valid for Caltech reconstruction) 
                % ======================================================================
                PLUME.T_plume                 = 1450;       % °C               
                PLUME.colat_lon_plume_center  = colat_lon_plume_center; % center of the plume (colatitude,longitude) [degrees]
                PLUME.r_plume                 = 100;        % radius of the plume [km]
                PLUME.l_TBC                   = 420;        % length of the plume [km] where Plume Temp BCs are imposed
                                                            % 420 km means 50 km inside the high-res region 
                                                            % (if depth high-res region = 300km)
                PLUME.d_TBC                   = 670;        % plume bottom depth [km] for imposing Temp Bcs
                PLUME.l_UBC                   = 420;        % length of the plume [km] where Plume Vel BCs are imposed
                PLUME.d_UBC                   = 670;        % plume bottom depth [km] for imposing Vel Bcs
                u_mean                        = plume_flux/(pi*PLUME.r_plume^2); % km/yr
                u_mean_scaled                 = u_mean * 1e6; % km/Myr
                PLUME.u_max                   = 2 * u_mean_scaled;
                PHYSICS.PLUME                 = PLUME;
            case 'no'
                PHYSICS.PLUME = [];
            otherwise
                error('SETTINGS.plume can olny be ''yes'' or ''no''')
        end
        
    case 'Debug'
%         % READ CRATONS INFO
%         % =================
%         pdir = pwd;
%         cd SETUP_TEST/GPLATES/CRATONS;
%         %------------------------------------------------------------------
%         % Check Cratons file format and number of files
%         %------------------------------------------------------------------
%         Cratons_info = dir();                % get info about the files in this directory
%         sample_file  = Cratons_info(4).name; % select 1st Vel_BCs file to check the extension
%         if strcmp(sample_file(end-2:end),'.xy')
%             system('./xy2txt.sh'); % change .xy extension to .txt extension
%         elseif strcmp(sample_file(end-3:end),'.txt')
%             % Cratons files are already in .txt format
%         else
%             error('wrong extension for CRATONS files. It should be .xy or .txt')
%         end
%         Cratons_info = dir('*.txt');
%         num_files    = numel(Cratons_info); % number of Vel_BCs files
%         
%         %------------------------------------------------------------------
%         % Create structure for CRATONS and load the files
%         %------------------------------------------------------------------
%         CRATONS(1:num_files) = struct('latlon',0);
%         fprintf('\n Loading Cratons info...');
%         for i = 1:num_files
%             fid               = fopen(Cratons_info(i).name,'r'); % open input file
%             tmp               = fscanf(fid, '%f');
%             CRATONS(i).latlon = reshape(tmp,2,[])';
%             fclose(fid); % close input file
%         end
%         fprintf(' done.\n\n');
%         cd(pdir);
%         PHYSICS.CRATONS = CRATONS;
%         clear Cratons_info i numfiles pdir sample_file
%     otherwise
%         error('typo in SETTINGS.plate_model')
end
%==========================================================================


%==========================================================================
% STOKES FLOW SOLUTION ALGORITHM
SETTINGS.method_eval_dens = 'interp_nodal';
  % Define how the density at each element's integration points (IPs) is
  % evaluated:
  % 'interp_nodal' --> will interpolate the nodal densities VAR.Dens at
  %                    integration points (without changing these
  %                    densities)
  % 'elem_phases'  --> will use the densities defined in PHYSICS.Dens and
  %                    element PhaseIDs: Dens = PHYSICS.Dens(MESH.PhaseID)
  % 'custom01'     --> will use the method defined for case 'custom_1' in 
  %                    the nested function "calc_dens_blk" (see flow
  %                    solver)
  % etc 
  % You can add more cases here but you have to edit the subfunction 
  % "calc_element_dens" (called during element assembly in the flow solver)

SETTINGS.method_eval_visc = 'interp_nodal';
  % Define how the viscosity at each element's integration points (IPs) 
  % is evaluated:
  % 'interp_nodal' --> will interpolate the nodal densities VAR.Visc at
  %                    integration points (without changing these
  %                    viscosities)
  % 'mean_nodal'   --> will calculate the harmonic mean of the nodal
  %                    viscosities for each element
  % 'elem_phases'  --> will use the viscosities defined in PHYSICS.Visc and
  %                    element PhaseIDs: Visc = PHYSICS.Visc(MESH.PhaseID)
  % 'custom01'     --> will use function "viscosity_custom01" to calculate
  %                    viscosity (dislocation creep only)
  % etc 
  % You can add more cases here but you have to edit the subfunction 
  % "calc_element_visc" (called during element assembly in the flow solver)
%==========================================================================


%==========================================================================
% CONSTANTS
PHYSICS.g    = -9.81;
  % [m/s^2] gravitational acceleration = -9.81; negative for downward !!
SETTINGS.gravity = 'radial';
  % Define the gravity direction
  % 'radial' --> radial direction (for spherical meshes)
  % 'z'      --> z direction (for cubic test meshes)
PHYSICS.Tm   = 1300;
  % [degree C] Temperature of mantle (initial mantle temperature)
PHYSICS.Ttop = 0; % 0;
  % [degree C] Temperature at top
PHYSICS.Tbot = 1300; % 0;
  % [degree C] Temperature at bottom
%==========================================================================


%==========================================================================
% MATERIAL PROPERTIES
PHYSICS.alpha    = 3.0e-5;
  % Thermal expansion coefficient (use 2.5e-5 to 3e-5)

PHYSICS.Dens     = 3300;
  % [kg/m^3] Reference density of mantle 
PHYSICS.minDens  = [];
PHYSICS.maxDens  = [];
  % [kg/m^3] Upper and lower limits for density (applies to all densities
  % calculated at integration points). Set empty [] to disable.

PHYSICS.Visc     = 2e18 / NUMSCALE.Visc0;
  % Reference viscosity of mantle [in units of NUMSCALE.Visc0]
PHYSICS.minVisc  = 1e18 / NUMSCALE.Visc0; % [in units of NUMSCALE.Visc0]
PHYSICS.maxVisc  = 1e23 / NUMSCALE.Visc0; % [in units of NUMSCALE.Visc0]
  % Upper and lower limits for viscosities (applies to all viscosities
  % calculated at integration points). Set empty [] to disable.

% Temperature dependence of mantle viscosity
PHYSICS.T_ref    = PHYSICS.Tm;
  % [degree C]; Reference temperature at which density is equal to 
  % PHYSICS.Dens and viscosity is equal to PHYSICS.Visc
PHYSICS.Ea  = 400e3; % 400e3; % 400e3;
  % [J/mol] activation energy; 300e3-600e3 J
  % Note: setting =0 disables temperature dependence

% Depth dependence of mantle viscosity
PHYSICS.Va  = 4e-6; %4e-6; %4e-6;
  % [m^3/mol] activation volume; 1e-6 to 1e-5
  % Note: setting =0 disables depth dependence
PHYSICS.a_UM2LM = 2500; %20;
  % factor viscosity increase in the lower mantle
    
% Thermal conductivity and specific heat capacity
PHYSICS.K        = 3.3;  % [W m^-1 K^-1]; provide in SI units !
PHYSICS.Cp       = 1000; % [J kg^-1 K^-1]; provide in SI units !
PHYSICS.kappa_SI = PHYSICS.K/(PHYSICS.Dens*PHYSICS.Cp); % diffusivity in SI
  % (using K=3.3 and Cp=1000 gives kappa=1e-6)
PHYSICS.kappa    = PHYSICS.kappa_SI / NUMSCALE.Kappa0; % scaled diffusivity
%==========================================================================


%==========================================================================
% THERMAL DIFFUSION
SETTINGS.T_solver = 'Diffusion';
  % Define method for solving diffusion thermal problem
  % 'Diffusion' --> Finite Element thermal solver for both 'std' and 'opt'
  %                 cases
  % 'Thermal'  --> Finite Element thermal solver (only 'opt' case: MILAMIN)
SETTINGS.OPTS_T.method = 'opt';
SETTINGS.OPTS_D.method = 'opt';
    % std = standard element assembly (slow)
    % opt = optimized assembly (blocks of elements; MILAMIN method)
    %       see Dabrowski et al. 2008 (doi:10.1029/2007GC001719)
SETTINGS.OPTS_T.method_eval_dens = SETTINGS.method_eval_dens; % use same as for Stokes flow
SETTINGS.OPTS_T.method_eval_cond = 'elem_phases';
SETTINGS.OPTS_T.method_eval_Cp   = 'elem_phases';
SETTINGS.OPTS_T.method_eval_dQdt = 'zero_dQdt';
  % Define how conductivity, heat capacity, and the energy source term 
  % are evaluated at each element's integration points
  % NOTE that "SETTINGS.method_eval_dens" (stokes flow setup above) is also
  % used in the thermal solver.
  % 'interp_nodal' --> will interpolate the nodal values (e.g. VAR.Cp) at
  %                    integration points (without changing these
  %                    densities)
  % 'elem_phases'  --> will use the values defined in, e.g., PHYSICS.Cp and
  %                    element PhaseIDs: Dens = PHYSICS.Cp(MESH.PhaseID)
  % 'elem_var'     --> will use element densities, e.g. VAR.Cp(iel)
  % method_eval_dQdt only:
  % 'zero_dQdt'    --> disables source term in energy equation
%==========================================================================

end % END OF FUNCTION model_parameters_sph