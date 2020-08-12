function M3TET_SPH_SOUTH_ATLANTIC_PLUME(plume_flux,colat_lon_plume_center)
% M3TET_SPH - 3D FINITE ELEMENT CONVECTION CODE IN SPHERICAL GEOMETRY
%
% Version M3TET_MELT also solves for melting rates in a multi-component
% (e.g peridotite matrix with embedded fertile or enriched veins) mantle.
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
% 
% M3TET solves for viscous flow and thermal advection/diffusion.
% Setup (physical parameters, etc) is defined in "model_parameters.m".
% Boundary conditions are defined in "bc_sph.m".
% Output is saved as mat-files in folder "SETTINGS.outdir".
% Tecplot 360 data files can be generated during post processing using
% function "BUILD_TEC360FILE_3D" located in POSTPROC.
%
% Patera's algorithm for solving the coupled pressure-velocity equations of
% Stokes flow. This algorithm uses a conjugate gradient formulation for 
% the pressure equation that is obtained after 'formally' inverting the
% momentum equations for velocity, and then substituting this into the
% incompressibilty equations. The velocity subproblem is solved using a
% Conjugate Gradient algorithm with a Cholesky direct solver on the
% coarsest mesh. A single V-cycle of geometric multigrid is used for
% preconditioning. M3TET uses the P2P1 Taylor-Hood element with continuous
% quadratic-order 10-node velocity and continuous linear-order
% 6-vertex-node pressure shape functions.
%
% Thermal advection diffusion is solved for by operator splitting:
% Semi-Lagrange advection scheme (Euler, 2nd-order Predictor-Corrector or 
% 4th-order Runge-Kutta back-tracking; linear, quadratic or cubic
% interpolation at back-tracking point).
% Finite element thermal diffusion solver (linear 4-node or quadratic
% 10-node elements; explicit/Crank-Nicholson/Galerkin/implicit time scheme;
% Conjugate Gradient solver with different preconditioners)
%
% JH/JPM Feb 2011
% JH Mar 2013: cleaned up and style made similar to POR3FLOW_SP
% JH Mar 2015: updated to latest M3TET version (e.g. added MUTILS)
% JH Feb 2016: updated to latest M3TET version
% JMT Jul 2016: Double Jacobian version (curved edges)
% 

clock0 = clock; t0 = tic; % initialize timing for profiling

fprintf(1,'\n\n\n\n');
fprintf(1,'=========================================================================\n');
fprintf(1,'            NEW CALCULATION USING M3TET_SPH_SOUTH_ATLANTIC\n');
fprintf(1,'=========================================================================\n');
addpath([pwd '/mfiles_M3TET']); % M3TET core functions
addpath([pwd '/mfiles_SPH/']);  % special functions for spherical coords
addpath([pwd '/SETUP_TEST/']);  % "model_parameter" - files are located here!
addpath([pwd '/TOOLS/subroutines_for_tests']); % functions for debugging
addpaths_mutils(); % add path to MUTILS installation on current computer

% SET THE DEFAULT ROOT RENDERER TO EITHER zbuffer OR OpenGL
% zbuffer makes nicer plots than OpenGL, but does not  support transparency
% zbuffer becomes very slow for larger 3D plots (e.g. a 3D scatter plot)
set(0, 'DefaultFigureRenderer', 'opengl');

%==========================================================================
% CHECK IF CODE RUNS IN PARALLEL
%==========================================================================
% The number of subdomains COMM.nsd is equal to number of Matlab "workers"
% or "labs" (command "numlabs"). The each subdomain's index COMM.myid is
% the index of each Matlab worker (command "labindex").
% If not running in parallel: COMM.nsd=1 and COMM.myid=1
COMM = init_COMM();

%==========================================================================
% DEFINING MODEL PARAMETERS
%==========================================================================
[PHYSICS,SETTINGS,NUMSCALE] = model_parameters_sph_south_atlantic_plume(plume_flux,colat_lon_plume_center);
ouput_on_screen = 0;

%==========================================================================
% CREATE OUTPUT FOLDER
%==========================================================================
if ~SETTINGS.restart
    create_output_folder(SETTINGS);
end

%======================================================================
% PREPARE LOGFILE OR OUTPUT IN TERMINAL
%======================================================================
SETTINGS = open_logfile(SETTINGS,COMM,ouput_on_screen);
fidl     = SETTINGS.fid_log;


if SETTINGS.restart
    %======================================================================
    % LOAD RESTART DATA AND CONTINUE AN EXISTING CALCULATION
    %======================================================================
    outdir = SETTINGS.outdir;
    [MESH,COMM,SETTINGS,PHYSICS,NUMSCALE,VAR,UBC,TBC,istep,iplot,time] = ...
        load_restart_data(SETTINGS,PHYSICS,COMM);
    SETTINGS.outdir    = outdir; clear outdir
    SETTINGS.fid_log   = fidl;
    SETTINGS.show_Figs = 'yes';

    % SAVE CONFIGURATION
    save([SETTINGS.outdir '/' COMM.prefix '_SETTINGS_restart'],...
         'SETTINGS','PHYSICS','NUMSCALE');

    fprintf(fidl,' Continuing calculation at time step %6i (time=%8.3f Myr) \n',...
        istep,time);
else
    %======================================================================
    % START A NEW CALCULATION
    %======================================================================

    % GLOBAL SWITCH TO DISABLE 'no' OR ENABLE 'yes' ALL FIGURES
    SETTINGS.show_Figs = 'yes';

    %======================================================================
    % LOAD 3D MESH FROM FILE, CREATE SUBDOMAINS, CRAETE MULTIGRID MESHES, 
    % AND SET UP INTER-SUBDOMAIN COMMUNICATION ARRAYS
    %======================================================================
    switch SETTINGS.edge_type
        case 'straight'
            [MESH,COMM] = init_mesh_straight_p(SETTINGS,PHYSICS,NUMSCALE,COMM);
        case 'curved'
            [MESH,COMM] = init_mesh_curved_p(SETTINGS,PHYSICS,NUMSCALE,COMM);
        otherwise
            error('SETTINGS.edge_type must be "straight" or "curved"')
    end
    
    % SAVE MESH AND COMM
    save([SETTINGS.outdir '/' COMM.prefix '_Mesh'],'MESH','-v7.3');
    save([SETTINGS.outdir '/' COMM.prefix '_COMM'],'COMM');
    
    %======================================================================
    % DISPLAY MODEL SIZE AND MESH PROPERTIES
    %======================================================================
    display_mesh_prop(MESH,COMM,SETTINGS,NUMSCALE);
    
    %======================================================================
    % INITIALIZATION OF VARIABLES
    %======================================================================
    [VAR,NUMSCALE] = init_fields_sph(MESH,SETTINGS,PHYSICS,NUMSCALE);
    
%     compare_serial_parallel_3d(MESH,COMM,VAR.T,'T',0);
%     compare_serial_parallel_3d(MESH,COMM,VAR.Dpl(:,1),'Dpl1',0);
%     compare_serial_parallel_3d(MESH,COMM,VAR.Dpl(:,2),'Dpl2',0);
%     compare_serial_parallel_3d(MESH,COMM,VAR.Vol(:,1),'Vol1',0);
%     compare_serial_parallel_3d(MESH,COMM,VAR.Vol(:,2),'Vol2',1);
    
    % SAVE CONFIGURATION
    save([SETTINGS.outdir '/' COMM.prefix '_SETTINGS'],...
         'SETTINGS','PHYSICS','NUMSCALE');
    
    %======================================================================
    % BOUNDARY CONDITIONS
    %======================================================================
    [UBC,TBC,VAR] = bc_sph(MESH,COMM,SETTINGS,PHYSICS,VAR,NUMSCALE);
    % SAVE BOUNDARY CONDITIONS
    save([SETTINGS.outdir '/' COMM.prefix '_BCs'],'UBC','TBC');
    
    if strcmp(SETTINGS.cratons_contours,'8') || strcmp(SETTINGS.cratons_contours,'8M')
        % Smooth shape of cratons -> Diffusion solver during 10 Myr with dt_diff = 0.1 Myr
        %----------------------------------------------------------------------------------
        dt_diff = 0.1; % time step for diffusion
        for i = 1:100
            VAR.T = diffusion3d_p(VAR.T,PHYSICS.kappa,MESH,COMM,SETTINGS,TBC,dt_diff);
        end
    end
end

%==========================================================================
% FEW MORE INITIALIZATIONS
%==========================================================================
%     switch SETTINGS.plate_model
%         case 'WJM'
%             varnames = {'Ux' 'Uy' 'Uz' 'Ue' 'Un' 'Ur' 'P' 'T' 'Dens' 'Visc'};
%         case 'GPlates'
%             varnames = {'Ux' 'Uy' 'Uz' 'Uth' 'Uph' 'Ur' 'P' 'T' 'Dens' 'Visc' 'Plate'};
%         case 'Debug'
%             varnames = {'Ux' 'Uy' 'Uz' 'Uth' 'Uph' 'Ur' 'P' 'T' 'Dens' 'Visc'};
%     end
if ~SETTINGS.restart
    time   = 0;
    istep  = 1;
    iplot  = 0;
    dt     = 0;
    time4plot = 0; % to create output in first time step
    
    % SAVE INITIAL FIELDS
    fprintf(fidl,' Writing output of initial fields\n');
    write_output_file(SETTINGS,VAR,COMM,-1,0,0,iplot);
%     write_tecplot_data_3d(MESH,VAR,SETTINGS,NUMSCALE,-1,0,0,varnames);
else
    % CALCULATE TIME STEP
    dt = calculate_dt(MESH,COMM,SETTINGS,NUMSCALE,VAR);
    
    % TIME FOR NEXT OUTPUT FILE
    time4plot = ceil((time-dt)/SETTINGS.save_period) * SETTINGS.save_period;
%     time4plot = (iplot-1) * SETTINGS.save_period;
end

% Initialize workspace for Semi-Lagrange advection solver
WS_SLM = [];

t_restart_period = 120; % write restart file every xxx sec
irestart         = 0;   % used to alternate between file #1 and #2
t_restart        = tic; % start timer for next restart file

fprintf(fidl,' INITIALIZATION      : %7.2f sec\n',toc(t0));

%==========================================================================
% TIME STEP LOOP
%==========================================================================
while time<SETTINGS.endtime && istep<=SETTINGS.nstep
    % DISPLAY TIME STEP, TIME, RUNTIME
    display_progress(istep,time,SETTINGS,NUMSCALE,clock0);
    
    % UPDATE VISCOSITY AND DENSITY OF MANTLE 
    VAR = rock_prop(VAR,PHYSICS,NUMSCALE,MESH,SETTINGS);

    % VARIABLE RANGES
    if mod(istep,1)==0 || istep==1
        display_variable_range_p(VAR,COMM,fidl);
    end

% 	% Add viscosity jump between upper and lower mantle
%     r           = sqrt(sum(MESH.GCOORD.^2,1));
%     VAR.Visc(:) = 1e21 / NUMSCALE.Visc0;
%     ind         = r>6370-670;
% %     VAR.Visc(ind) = 5e21 / NUMSCALE.Visc0;
% %     ind           = r>410;
%     VAR.Visc(ind) = 1e19 / NUMSCALE.Visc0;

%     % Add density anomaly to drive internal flow
%     VAR.Dens(MESH.GCOORD(3,:)>5000) = 3310;
    
    % UPDATE VELOCITY BOUNDARY CONDITIONS (only for SETTINGS.plate_model = 'GPlates')
    if strcmp(SETTINGS.plate_model,'GPlates')
        [UBC,VAR] = update_vel_bc(MESH,COMM,PHYSICS,VAR,UBC,time);
    end
    % SOLVE FOR VISCOUS FLOW (VELOCITY AND PRESSURE)
    VAR = flow3d_patera_sph(VAR,MESH,COMM,UBC,SETTINGS,PHYSICS,NUMSCALE,dt);

%     if numlabs==1
%         % Integrate flow across core-mantle boundary...
%         nods_bot = find(MESH.PointID==301);
%         NN_INTEG = boundary_integration_matrix_3d(MESH.GCOORD,MESH.EL2NOD{1},nods_bot);
%         flux_bot = sum(NN_INTEG * VAR.Ur);
%         area_bot = sum(NN_INTEG * ones(MESH.nnod,1));
% 
%         % ... and top surface
%         nods_top = find(MESH.PointID==306);
%         NN_INTEG = boundary_integration_matrix_3d(MESH.GCOORD,MESH.EL2NOD{1},nods_top);
%         flux_top = sum(NN_INTEG * VAR.Ur);
%         area_top = sum(NN_INTEG * ones(MESH.nnod,1));
%         
%         fprintf(' Fluxes across CMB:\n area / total flux / flux per area\n %.4e / %.4e / %.4e\n',...
%             area_bot,flux_bot,flux_bot/area_bot);
%         fprintf(' Fluxes across surface:\n area / total flux / flux per area\n %.4e / %.4e / %.4e\n',...
%             area_top,flux_top,flux_top/area_top);
%     end
    
	% CALCULATE TIME STEP
    dt = calculate_dt(MESH,COMM,SETTINGS,NUMSCALE,VAR);

    % WRITE MATLAB DATA FILE
    if (time4plot-time)<dt || time+dt>=SETTINGS.endtime
        iplot     = iplot + 1;
%         write_tecplot_data_3d(MESH,VAR,SETTINGS,NUMSCALE,time,istep,iplot,varnames);
        write_output_file(SETTINGS,VAR,COMM,time,istep,dt,iplot);
        time4plot = time4plot + SETTINGS.save_period;
%         time4plot = iplot * SETTINGS.save_period;
    end
    
    speed     = sqrt(VAR.Ux.^2+VAR.Uy.^2+VAR.Uz.^2);
    rms_speed = rms(speed);
    rms_speed
    max_speed = max(speed);
    max_speed
    
%     compare_serial_parallel_3d(MESH,COMM,VAR.T,'T0',0);
    
    % SOLVE FOR THERMAL ADVECTION-DIFFUSION
    [VAR,WS_SLM] = thermal_advdiff_3d_p...
        (VAR,MESH,COMM,PHYSICS,SETTINGS,NUMSCALE,TBC,dt,WS_SLM);
    
    % UPDATE TIME AND TIME STEP COUNTER
    time  = time + dt;
    istep = istep + 1;
    
    % WRITE RESTART FILE
    if toc(t_restart)>t_restart_period
        irestart  = mod(irestart,2)+1;
        write_restart_file(SETTINGS,VAR,COMM,istep,iplot,time,dt,irestart);
        t_restart = tic;
    end
end
%==========================================================================
% END OF TIME STEP LOOP
%==========================================================================

% DISPLAY TIME STEP, TIME, RUNTIME
display_progress(istep,time,SETTINGS,NUMSCALE,clock0);

fprintf(1, '\n');
fprintf(1, ' CALCULATION FINISHED. Output files in folder:\n "%s"\n',...
    SETTINGS.outdir);
fprintf(1,'=========================================================================\n');

end % END OF FUNCTION M3TET_SPH