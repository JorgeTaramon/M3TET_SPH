function [MESH,COMM,SETTINGS,PHYSICS,NUMSCALE,VAR,UBC,TBC,istep,iplot,time] = ...
              load_restart_data(SETTINGS,PHYSICS,COMM)

output_dir = [SETTINGS.outdir '/'];
fprintf(' Loading data from folder\n   "%s\n to restart the calculation\n',...
        output_dir);
% Check if output folder exists
a = ls(output_dir);
if isempty(a)
    error('Folder "%s" does not exist! Cannot restart.',output_dir);
end

if isfield(SETTINGS,'nsd_restart') && SETTINGS.nsd_restart>0 && ...
    SETTINGS.nsd_restart~=COMM.nsd
   [MESH,VAR,COMM,SETTINGS,PHYSICS,time,istep,iplot] = ...
       load_data_different_nsd(output_dir,COMM,SETTINGS,PHYSICS);
else
   [MESH,VAR,UBC,TBC,COMM,SETTINGS,PHYSICS,NUMSCALE,time,istep,iplot] = ...
       load_data_same_nsd(output_dir,COMM,SETTINGS,PHYSICS);
end

% Initialize COMM again (a calculation may have been started on a unix
% machine and is now continued on windows, or vice versa; in this case the
% number of threads for MUTILS will change)
COMM = init_COMM();
 
end % END OF FUNCTION load_restart_data

% #########################################################################
%                               SUBFUNCTIONS
% #########################################################################


function [MESH,VAR,PBC,TBC,COMM,SETTINGS,PHYSICS,time,istep,iplot] = ....
             load_data_different_nsd(output_dir,COMM,SETTINGS,PHYSICS)
error(' Not yet working (boundary conditions...');

% Load "MESH"
fprintf(' Loading "MESH"\n ');
nsd_restart = SETTINGS.nsd_restart;
GCOORD_rst  = cell(1,nsd_restart);
ECOORD_rst  = cell(1,nsd_restart);
prefix_rst  = cell(1,nsd_restart);
for isd=1:nsd_restart
    prefix_rst{isd} = ['Lab' num2str_d(nsd_restart,2) 'x' num2str_d(isd,2)];
    file            = [prefix_rst{isd} '_Mesh.mat'];
    fprintf('from file %s...\n ',file);
    LOADED          = load([output_dir file],'MESH');
    GCOORD_rst{isd} = single(LOADED.MESH.GCOORD);
    ECOORD_rst{isd} = single(calc_tetra_center(LOADED.MESH.GCOORD,LOADED.MESH.EL2NOD{1}));
end
fprintf('done\n');

%======================================================================
% LOAD 3D MESH FROM FILE, CREATE SUBDOMAINS, AND
% SET UP INTER-SUBDOMAIN COMMUNICATION ARRAYS
%======================================================================
[MESH,COMM] = init_mesh_p(SETTINGS,PHYSICS,COMM);
% Save mesh data in output folder
save([output_dir COMM.prefix '_Mesh'],'MESH');
save([output_dir COMM.prefix '_Comm'],'COMM');
labBarrier

fprintf(' Calculating pointers between generated and loaded mesh...');
ptr_nod_D  = cell(1,nsd_restart);
ptr_nod_SD = cell(1,nsd_restart);
check_nod  = zeros(MESH.nnod,1,'int32');
ptr_elm_D  = cell(1,nsd_restart);
ptr_elm_SD = cell(1,nsd_restart);
check_elm  = zeros(MESH.nnod,1,'int32');
ECOORD     = single(calc_tetra_center(MESH.GCOORD,MESH.EL2NOD{1}));
for isd=1:nsd_restart
    [~,ia,ib] = intersect(single(MESH.GCOORD'),GCOORD_rst{isd}','rows');
    ptr_nod_D {isd} = ia; % VAR.T(ia) = DATA_SD.VAR.T(ib);
    ptr_nod_SD{isd} = ib; %
    check_nod(ia)   = 1;
    [~,ia,ib] = intersect(ECOORD',ECOORD_rst{isd}','rows');
    ptr_elm_D {isd} = ia; % VAR.T(ia) = DATA_SD.VAR.T(ib);
    ptr_elm_SD{isd} = ib; %
    check_elm(ia)   = 1;
end
if any(~check_nod)
    error(' Could not find all nodes.');
end
if any(~check_elm)
    error(' Could not find all elements.');
end
fprintf('done\n');

if strcmp(SETTINGS.restart_setup,'from_file')
    % Load "SETTINGS" and "PHYSICS"
    fprintf(' Loading "PHYSICS" and "SETTINGS"...');
    try
        RestartData = load([output_dir prefix_rst{1} '_SETTINGS.mat']);
    catch
        RestartData = load([output_dir 'SETTINGS.mat']);
    end
    PHYSICS                = RestartData.PHYSICS;
    Settings0              = SETTINGS;
    SETTINGS               = RestartData.SETTINGS;
    SETTINGS.restart       = 1;
    SETTINGS.restart_file  = Settings0.restart_file;
    SETTINGS.restart_setup = Settings0.restart_setup;
    SETTINGS.time_end      = Settings0.time_end;
    SETTINGS.outdir        = Settings0.outdir;
    fprintf('done\n');
end

% SAVE CONFIGURATION
save([SETTINGS.outdir '/' COMM.prefix '_SETTINGS'],'SETTINGS','PHYSICS');

% Load "VAR"
fprintf(' Loading "VAR"\n ');
VAR_rst = cell(1,nsd_restart);
if isnumeric(SETTINGS.restart_file) || str2double(SETTINGS.restart_file)>-1
    if ischar(SETTINGS.restart_file)
        SETTINGS.restart_file = str2double(SETTINGS.restart_file);
    end
    try
        for isd=1:nsd_restart
            file = [prefix_rst{isd} '_Data_' num2str_d(SETTINGS.restart_file,4) '.mat'];
            [VAR_rst{isd},time,istep] = load_restart_file([output_dir '/' file]);
            fprintf('from file %s...\n',file);
        end
    catch %#ok<CTCH>
        error('Could not open file "%s" in folder "%s"',...
              file,output_dir);
    end
    iplot = 1 + SETTINGS.restart_file;
    
elseif ischar(SETTINGS.restart_file) && strcmp(SETTINGS.restart_file,'last')
    RestartFiles = ls([output_dir prefix_rst{1} '_RestartFile*.mat']);
    if isempty(RestartFiles)
        error('No restart files in output folder!');
    end
    RestartFiles = dir([output_dir prefix_rst{1} '_RestartFile*.mat']);
    if length(RestartFiles)==1
        file = RestartFiles(1).name;
    else
        if RestartFiles(1).datenum>RestartFiles(2).datenum
            file = RestartFiles(1).name;
        else
            file = RestartFiles(2).name;
        end
    end
    for isd=1:nsd_restart
        file(1:8) = prefix_rst{isd};
        fprintf('from file %s...\n ',file);
        [VAR_rst{isd},time,istep,iplot] = load_restart_file([output_dir '/' file]);
    end
else
    error('Invalid value in "SETTINGS.restart_file"');
end
fprintf('done\n');

fprintf(' Assembling VAR from restart data...\n ');
for isd=1:nsd_restart
    varlist  = fieldnames(VAR_rst{isd});
    var_type = -ones(1,length(varlist));
    for ivar=1:length(varlist)
        varname = varlist{ivar};
        if length(VAR_rst{isd}.(varname))==size(GCOORD_rst{isd},2)
            var_type(ivar) = 1;
        elseif length(VAR_rst{isd}.(varname))==size(ECOORD_rst{isd},2)
            var_type(ivar) = 0;
        else
            error(' Unknown variable type.');
        end
        if var_type(ivar)
            if isd==1
                VAR.(varname) = zeros(MESH.nnod,1);
            end
            VAR.(varname)(ptr_nod_D{isd}) = VAR_rst{isd}.(varname)(ptr_nod_SD{isd});
        else
            if isd==1
                VAR.(varname) = zeros(MESH.nel,1);
            end
            VAR.(varname)(ptr_elm_D{isd}) = VAR_rst{isd}.(varname)(ptr_elm_SD{isd});
        end
    end
end
fprintf('done\n');

error('This block needs to be verified.');
% Load "BCs"
fprintf(' Loading "BCs"...');
for isd=1:nsd_restart
    file    = [prefix_rst{isd} '_BCs.mat'];
    fprintf('from file %s...\n ',file);
    LOADED  = load([output_dir file]);
    ind     = ismember(ptr_nod_SD{isd},LOADED.PBC.nod);
    PBC.nod = ptr_nod_D{isd}( LOADED.PBC.nod(ptr_nod_SD{isd}(ind)) );
    PBC.val = LOADED.PBC.val(ptr_nod_SD{isd}(ind));
    if isfield(LOADED.PBC,'Ptop')
        PBC.Ptop = LOADED.PBC.Ptop;
    end
    if isfield(LOADED.PBC,'Pbot')
        PBC.Pbot = LOADED.PBC.Pbot;
    end
end
TBC = bc_temperature(MESH,SETTINGS,PHYSICS,VAR);
fprintf('done\n');

end % END OF SUBFUNCTION load_data_different_nsd

% #########################################################################

function [MESH,VAR,UBC,TBC,COMM,SETTINGS,PHYSICS,NUMSCALE,time,istep,iplot] = ....
             load_data_same_nsd(output_dir,COMM,SETTINGS,PHYSICS)

% Load "MESH"
fprintf(' Loading "MESH"...');
RestartData = load([output_dir COMM.prefix '_Mesh.mat']);
MESH        = RestartData.MESH;
fprintf('done\n');

% % Load "COMM"
% fprintf(' Loading "Communication data"...');
% RestartData = load([output_dir COMM.prefix '_Comm.mat']);
% COMM        = RestartData.COMM;
% fprintf('done\n');

% Load "NUMSCALE"
fprintf(' Loading "NUMSCALE"...');
RestartData = load([output_dir COMM.prefix '_SETTINGS.mat']);
NUMSCALE    = RestartData.NUMSCALE;
fprintf('done\n');

if strcmp(SETTINGS.restart_setup,'from_file')
    % Load "SETTINGS" and "PHYSICS"
    fprintf(' Loading "PHYSICS" and "SETTINGS"...');
    try
        RestartData = load([output_dir COMM.prefix '_SETTINGS.mat']);
    catch
        RestartData = load([output_dir 'SETTINGS.mat']);
    end
    PHYSICS                = RestartData.PHYSICS;
    Settings0              = SETTINGS;
    SETTINGS               = RestartData.SETTINGS;
    SETTINGS.restart       = 1;
    SETTINGS.restart_file  = Settings0.restart_file;
    SETTINGS.restart_setup = Settings0.restart_setup;
    SETTINGS.endtime       = Settings0.endtime;
    fprintf('done\n');
end

% Load "BCs"
fprintf(' Loading "BCs"...');
RestartData = load([output_dir COMM.prefix '_BCs.mat']);
UBC         = RestartData.UBC;
TBC         = RestartData.TBC;
fprintf('done\n');

% Load "VAR"
fprintf(' Loading "VAR" ');
if isnumeric(SETTINGS.restart_file) || str2double(SETTINGS.restart_file)>-1
    if ischar(SETTINGS.restart_file)
        SETTINGS.restart_file = str2double(SETTINGS.restart_file);
    end
    try
        file        = [COMM.prefix '_Data_' num2str_d(SETTINGS.restart_file,4) '.mat'];
        fprintf('from file %s...',file);
        [VAR,time,istep] = load_restart_file([output_dir '/' file]);
    catch
        error('Could not open file "%s" in folder "%s"',...
              file,output_dir);
    end
    iplot = SETTINGS.restart_file;
    
elseif ischar(SETTINGS.restart_file) && strcmp(SETTINGS.restart_file,'last')
    RestartFiles = ls([output_dir COMM.prefix '_RestartFile*.mat']);
    if isempty(RestartFiles)
        error('No restart files in output folder!');
    end
    RestartFiles = dir([output_dir COMM.prefix '_RestartFile*.mat']);
    if length(RestartFiles)==1
        file = RestartFiles(1).name;
    else
        if RestartFiles(1).datenum>RestartFiles(2).datenum
            file = RestartFiles(1).name;
        else
            file = RestartFiles(2).name;
        end
    end
    fprintf('from file %s...',file);
    [VAR,time,istep,iplot] = load_restart_file([output_dir '/' file]);
else
    error('Invalid value in "SETTINGS.restart_file"');
end
fprintf('done\n');

end

% #########################################################################

function [VAR,time,istep,iplot] = load_restart_file(filename)

RestartData = load(filename);
if isfield(RestartData,'time_yr')
    time = RestartData.time_yr;
else
    time = RestartData.time;
end
istep = RestartData.istep;
VAR   = RestartData.VAR;
if isfield(RestartData,'iplot')
    iplot = RestartData.iplot;
else
    iplot = [];
end

end % END OF SUBFUNCTION load_restart_file