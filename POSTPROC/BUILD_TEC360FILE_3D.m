function BUILD_TEC360FILE_3D(OPTS)

pdir=pwd;cd('..');addpath([pwd '/mfiles_M3TET']);cd(pdir);

% DEFAULTS
% ========
% (1) specify output directory (folder where data is located)
% outdir = 'D:\Projects\SPHERE\Trash2';
outdir = '/Users/jorge/Tests/SPH/Trash_01';
% (2) specify which output files are converted
%     (number OR 'all' OR 'last')
% plotfiles = [0 113];
plotfiles = 'all';
% plotfiles = 'last';
% plotfiles = 74;
% (3) specify total number of subdomains
%     empty ([]) will find out looking in outdir
nsd = 1; %[];
% (4) specify which subdomain data is converted
%     empty ([]) will convert all subdomains
indSD = [];
% (5) specify which variables are converted
var_choice = 'all';
% var_choice = {'Vz' 'T'};
% (6) specify suffix of Tecplot file name
%     Tecplot file will be named ['Lab32x24' suffix]
suffix = '_Tec360.dat';
% ========

% Overwrite defaults if different values are specified in structure OPTS
if nargin>0
    if isfield(OPTS,'outdir')
        outdir = OPTS.outdir;
    end
    if isfield(OPTS,'plotfiles')
        plotfiles = OPTS.plotfiles;
    end
    if isfield(OPTS,'nsd')
        nsd = OPTS.nsd;
    end
    if isfield(OPTS,'indSD')
        indSD = OPTS.indSD;
    end
    if isfield(OPTS,'var_choice')
        var_choice = OPTS.var_choice;
    end
    if isfield(OPTS,'suffix')
        suffix = OPTS.suffix;
    end
end

if iscell(var_choice)
    vars_included = var_choice;
else
    switch var_choice  % X, Y, Z will always be selected !!
        % SETUPS FOR 3D HYDROTHERMAL CONVECTION DATA
        case 'HT_few'
            vars_included = {'T' 'P' 'Vx' 'Vy' 'Vz' 'Rho_f'};
        case 'HT_standard'
            vars_included = {'T' 'P' 'Vx' 'Vy' 'Vz' 'Perm' 'Rho_f' 'Cp_f' 'Mu_f'};
        % SETUPS FOR 3D MANTLE CONVECTION DATA
        case 'M3_few'
            vars_included = {'Ux' 'Uy' 'Uz' 'T'};
        case 'M3_standard'
            vars_included = {'Ux' 'Uy' 'Uz' 'T' 'P' 'Dens' 'Visc'};
        case 'M3_melting'
            vars_included = {'Ux' 'Uy' 'Uz' 'T' 'P' 'Dens' 'Visc' 'M' 'Dpl' 'Vol'};
        % SETUPS FOR ANY 3D DATA
        case 'onlyT'
            vars_included = {'T'};
        case 'all'
            vars_included = 'all';
        otherwise
            error(' Unknown variable choice');
    end
end

% -------------------------------------------------------------------------
%      No editing of this file should be required below this line
% -------------------------------------------------------------------------

outdir  = [outdir filesep];
workdir = define_working_directory(outdir); % *SUBFUNCTION*

if isempty(nsd)
    [indSD,nsd] = get_sd_list(outdir,indSD);
else
    if isempty(indSD)
        indSD = 1:nsd;
    end
end
if numlabs>1
    indSD_all = indSD;
    indx      = labindex:numlabs:length(indSD_all);
    indSD     = indSD_all(indx);
end

fprintf(' Will processing data of SDs:');
fprintf(' %2i',indSD);fprintf('\n');

for ii=1:length(indSD)
    isd     = indSD(ii);
    fprintf(' Processing data of SD %1i of %1i\n',isd,length(indSD));

    prefix  = ['Lab' num2str_d(nsd,2) 'x' num2str_d(isd,2)];
    tecfile = [prefix suffix];
    
    % --------------------------------------------------------------------
    matfiles = get_matfile_list(outdir,[prefix '_Data_'],plotfiles); % *SUBFUNCTION*
    nplot    = length(matfiles);
    cell_out = cell(1,4+2*(nplot-1));
    % --------------------------------------------------------------------

    fprintf(' Loading numerical scaling data...');
    try
        data    = load([outdir prefix '_SETTINGS.mat']);
        Pscale  = 1e-9 * data.NUMSCALE.P0;
        Muscale = data.NUMSCALE.Visc0;
        clear data
    catch %#ok<CTCH>
        Pscale  = 1;
        Muscale = 1;
    end
    fprintf(' done\n');

    % Load MESH data
    fprintf(' Loading mesh data...');
    try
        data  = load([outdir prefix '_Mesh.mat']);
    catch %#ok<CTCH>
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
    
    GCOORD = MESH.GCOORD;
    EL2NOD = MESH.EL2NOD{1};
    nnodel = size(EL2NOD,1);
    nnod   = max(max(EL2NOD));
    nVnod  = max(max(EL2NOD(1:4,:)));
    if nnodel==10
        EL2NOD_10 = EL2NOD;
        EL2NOD    = tetmesh_p2_to_p1(GCOORD,EL2NOD);
    end
    nel    = size(EL2NOD,2);
    fprintf(' done\n');
    
    icell = 1; new_file = 1;
    for ifile=1:nplot
        fprintf(' Processing data file #%1i of %1i\n',ifile,nplot);
        matfile = matfiles{ifile};
        try
            data = load(matfile);
        catch %#ok<CTCH>
            error(' Could no open file "%s"\n',matfile);
        end
        
        if isfield(data,'Var')
            VAR = data.Var;
        else
            VAR = data.VAR;
        end
        
        P_vertices     = VAR.P;
        VAR.P          = zeros(MESH.nnod,1);
        VAR.P(1:nVnod) = P_vertices;
        VAR.P(MESH.EL2NOD{1}(5:10,:)) = ...
            0.5*(P_vertices(MESH.EL2NOD{1}([1 2 3 4 1 2],:)) ...
               + P_vertices(MESH.EL2NOD{1}([2 3 4 1 3 4],:)));
        clear P_vertices
        VAR.P    = VAR.P * Pscale;
        VAR.Visc = log10(VAR.Visc * Muscale);
        
        if isfield(VAR,'Z')
            zcoords_change = 1;
        else
            zcoords_change = 0;
        end
        
        if isfield(data,'time')
            time = data.time;
        else
            time = data.time_yr;
        end
        clear data
        
        if ifile==1
            zonename  = sprintf('Time = %.3f',time);
            fprintf(' Will write the following variables to file:\n X Y Z ');
            var_names = get_variable_names(); % * NESTED FUNCTION*
            fprintf('\n');
            nvar      = length(var_names);
            [cell_out{icell},mvar] = make_tec360_varlist(); % * NESTED FUNCTION*
            icell                  = icell + 1;
        end

        cell_out{icell} = make_tec360_header(ifile); % * NESTED FUNCTION*
        icell           = icell + 1;

        % Write data
        cell_out{icell} = make_tec360_datablock(ifile); % * NESTED FUNCTION*
        icell           = icell + 1;

        if ifile==1
            write_tec360_ascii_file(); % * NESTED FUNCTION*
            clear cell_out
            icell = 1; new_file = 0;
            cell_out = cell(1,4+2*(nplot-1));
            
            % write connectivity
            cell_out{icell} = sprintf('%1i %1i %1i %1i\n',EL2NOD);
            icell           = icell + 1;
            clear EL2NOD
            
            write_tec360_ascii_file(); % * NESTED FUNCTION*
            clear cell_out
            icell = 1; new_file = 0;
            cell_out = cell(1,4+2*(nplot-1));
        end

        tmp = whos('cell_out'); MB = 0.5*tmp.bytes/1024^2;
        if MB>500 || ifile==nplot && icell>1
            write_tec360_ascii_file(); % * NESTED FUNCTION*
            clear cell_out
            icell = 1; new_file = 0;
            cell_out = cell(1,4+2*(nplot-1));
        end
    end

    fprintf(' Tecplot data file "%s" was written to folder\n "%s"\n',...
        tecfile,workdir);
    fprintf(' It contains %1i time steps and the following variables\n',nplot);
    if zcoords_change
        disp(['X' 'Y' 'Z' var_names]);
    else
        disp(var_names);
    end
    make_tec360_binary_file(); % * NESTED FUNCTION*
end

% =========================================================================
%                            NESTED FUNCTIONS
% =========================================================================

function write_tec360_ascii_file()
    tmp = whos('cell_out'); MB = 0.5*tmp.bytes/1024^2;
    fprintf(' Writing %5.1fMB to tecfile...',MB);
    if new_file
        output_unit = fopen([workdir '/' tecfile],'w');
    else
        output_unit = fopen([workdir '/' tecfile],'a');
    end
    % Write tecplot file
    if MB<1500
        fprintf(output_unit, [cell_out{1:icell-1}]);
    else
        for i=1:icell-1
            fprintf(output_unit, [cell_out{i}]);
        end
    end
    fclose(output_unit);
    fprintf('done.\n');
end % END OF NESTED FUNCTION write_tec360_ascii_file

% =========================================================================

function varnames = get_variable_names()
    all_vars = fieldnames(VAR);
    varnames = cell(0);
    for ivar=1:length(all_vars)
        select_var = 0;
        if length(VAR.(all_vars{ivar}))==nnod || ...
           length(VAR.(all_vars{ivar}))==nVnod
            if strcmp(vars_included,'all')
                select_var = 1;
            else
                for kk=1:length(vars_included)
                    if strcmp(all_vars{ivar},vars_included{kk})
                        select_var = 1;
                        break
                    end
                end
            end
        end
        if select_var
            varnames{end+1} = all_vars{ivar}; %#ok<AGROW>
            fprintf(' %s ',all_vars{ivar});
        end
    end
end % END OF NESTED FUNCTION get_variable_names

% =========================================================================

function [s,mvar] = make_tec360_varlist()
    s    = 'VARIABLES = "X", "Y", "Z", ';
    mvar = 0;
    for ivar=1:nvar
        var_name = var_names{ivar};
        m = size(VAR.(var_name),2);
        if m==1
            s = [s sprintf('"%s", ',var_name)]; %#ok<AGROW>
        else
            for ic=1:m
                s = [s sprintf('"%s%1i", ',var_name,ic)]; %#ok<AGROW>
            end
        end
        mvar = mvar + m;
    end
    s(end-1:end) = '\n';
end % END OF NESTED FUNCTION make_tec360_varlist

% =========================================================================

function s = make_tec360_header(ifile)
    if ifile==1
        s = sprintf('ZONE T="%s",N=%1i, E=%1i, DataPacking=POINT, ZoneType=FETETRAHEDRON, SolutionTime=%8.4f \n',...
              zonename,nnod,nel,time);
    else
        if zcoords_change
            s = sprintf('ZONE T="%s",N=%1i, E=%1i, DataPacking=POINT, ZoneType=FETETRAHEDRON, VARSHARELIST = ([1, 2]=1), CONNECTIVITYSHAREZONE = 1, SolutionTime=%8.4f \n',...
                  zonename,nnod,nel,time);
        else
            s = sprintf('ZONE T="%s",N=%1i, E=%1i, DataPacking=POINT, ZoneType=FETETRAHEDRON, VARSHARELIST = ([1, 2, 3]=1), CONNECTIVITYSHAREZONE = 1, SolutionTime=%8.4f \n',...
                  zonename,nnod,nel,time);
        end
    end
end % END OF NESTED FUNCTION make_tec360_header

% =========================================================================

function s = make_tec360_datablock(ifile)
    if ifile==1
        varblock        = zeros(3+mvar,nnod); % xyz-coord
        varblock(1:3,:) = MESH.GCOORD;
        if zcoords_change && isfield(VAR,'Z')
            varblock(3,:) = VAR.Z(:)';
        end
        i1              = 4;
    else
        if zcoords_change
            varblock      = zeros(1+mvar,nnod); % z-coord
            varblock(1,:) = VAR.Z(1:nnod);
            i1            = 2;
        else
            varblock = zeros(mvar,nnod);
            i1       = 1;
        end
    end
    
    for ivar=1:nvar
        var_name = var_names{ivar};
        var      = VAR.(var_name)';
        
        if size(var,2)<nnod
            if size(var,2)==nVnod
                var = interp_edge_node_vals(EL2NOD_10,var);
            else
                error(' Dont understand size of variable "%s"',var_name);
            end
        end
        
        i2 = i1 + size(var,1) - 1;
        if size(var,1)~=i2-i1+1 || size(var,2)~=size(varblock,2)
            error('dimensions  do not fit');
        end
        varblock(i1:i2,:) = var;
        i1 = i2 + 1;
    end

    fmt = repmat('%6E ',1,size(varblock,1));
    fmt(end:end+1) = '\n';
    s   = sprintf(fmt, varblock);
end % END OF NESTED FUNCTION make_tec360_datablock

% =========================================================================

function make_tec360_binary_file()
    fprintf('\n Converting Tecplot360 file from ASCII to binary...');
    pdir=pwd;cd(workdir);
    tecfile2 = [tecfile(1:end-3) 'plt'];
    [~,result] = system(['preplot ' tecfile ' ' tecfile2]);
    if isempty(strfind(result,'ERROR')) && ...
       isempty(strfind(result,'error')) && ...
       isempty(strfind(result,'Err: Cannot find'))
        fprintf('done.\n Binary file named "%s"\n',tecfile2);
        fprintf(' Deleting ascii file...');
        cmd = sprintf('del /Q /F "%s"',tecfile);
        [status,result] = system(cmd);
        if status || ~isempty(strfind(result,'cannot'))
            disp(result);
        else
            fprintf('done.\n');
        end
    else
        fprintf('FAILED.\n %s\n',result);
    end

    if ~isempty(workdir)
        fprintf(' Moving binary file(s) to data directory\n "%s"\n...',outdir);
        cmd = sprintf('move /Y "%s\\%s" "%s\\"',workdir,tecfile2,outdir(1:end-1));
        [status,result] = system(cmd);
        if status
            disp(result);
        else
            fprintf('done.\n');
        end
    end
    cd(pdir);
end % END OF NESTED FUNCTION make_tec360_binary_file

end % END OF FUNCTION BUILD_TEC360FILE_3D

% #########################################################################
%                               SUBFUNCTIONS
% #########################################################################

function workdir = define_working_directory(outdir)

tmp      = pctconfig;
hostname = tmp.hostname;
switch hostname
    case 'b3pc16'
        workdir = 'C:\Users\jhasenclever\scratch-SSD';
    case 'b3pc17'
        workdir = 'E:\PLT-conversion';
    case 'Maupiti2'
        workdir = 'C:\Users\joha\scratch-SSD';
    otherwise
        fprintf(' No fast (SSD) working direktory has been defined for this computer.\n');
        fprintf(' If this computer has a SSD, edit subfunction "define_working_directory".\n');
        workdir = '';
end

pdir = pwd;
if isempty(workdir)
    workdir = outdir;
else
    workdir = [workdir '/Lab' num2str_d(numlabs,2) 'x' num2str_d(labindex,2) '/'];
    try
        cd(workdir);
    catch %#ok<CTCH>
        [~,~,~] = mkdir(workdir);
    end
end
cd(pdir);

end % END OF SUBFUNCTION define_working_directory

% #########################################################################

function [indSD,nsd] = get_sd_list(outdir,indSD)

content = dir(outdir);
if isempty(content)
    error('Folder "%s" does not exist.',outdir);
end

nsd = [];
for j=1:length(content)
    filename = content(j).name;
    if length(filename)>11 && ...
       strcmp(filename(1:3),'Lab') &&...
       strcmp(filename(end-8:end),'_MESH.mat')
       nsd(end+1) = str2double(filename(4:5)); %#ok<AGROW>
    end
end
nsd = unique(nsd);
if length(nsd)~=1
    error(' Cannot find out how many CPUs produced output in folder\n"%s"',outdir');
end
if isempty(indSD)
    indSD = 1:nsd;
else
    ind = ~ismember(indSD,1:nsd);
    if any(ind)
        error(' This subdomain does not exist in data folder: %1i\n',indSD(ind));
    end
end

end % END OF SUBFUNCTION get_sd_list

% #########################################################################

function matfiles = get_matfile_list(outdir,datafile_pattern,plotfiles)

content = dir(outdir);
if isempty(content)
    error('Folder "%s" does not exist.',outdir);
end

alldatafiles = cell(0);
ind_plot     = [];
n = length(datafile_pattern);
for i=1:length(content)
    filename = content(i).name;
    if length(filename)>n && ...
       strcmp(filename(1:n),datafile_pattern) && ...
       strcmp(filename(end-3:end),'.mat')
        alldatafiles{end+1} = [outdir filename]; %#ok<AGROW>
        j               = strfind(filename,'.mat')-1;
        ind_plot(end+1) = str2double(filename(n+1:j)); %#ok<AGROW>
    end
end
if isempty(ind_plot)
    error('Cannot find Matlab data files');
end
if ischar(plotfiles)
    switch plotfiles
        case 'last'
            matfiles = alldatafiles(end);
        case 'all'
            matfiles = alldatafiles;
        otherwise
            error('Invalid value of plotfiles');
    end
else
    matfiles = cell(1,length(plotfiles)); j = 0;
    for i=1:length(ind_plot);
        if ismember(ind_plot(i),plotfiles)
            j           = j + 1;
            matfiles{j} = alldatafiles{i};
        end
    end
    matfiles(j+1:end) = [];
end

end % END OF SUBFUNCTION get_matfile_list

% #########################################################################

function var = interp_edge_node_vals(EL2NOD,var0)

nnod  = max(EL2NOD(:));
nVnod = size(var0,2);
nc    = size(var0,1);
var   = zeros(nc,nnod);
var(:,1:nVnod) = var0;
for ic=1:nc
    var(ic,EL2NOD( 5,:)) = 0.5.*(var(ic,EL2NOD(1,:)) + var(ic,EL2NOD(2,:)));
    var(ic,EL2NOD( 6,:)) = 0.5.*(var(ic,EL2NOD(2,:)) + var(ic,EL2NOD(3,:)));
    var(ic,EL2NOD( 7,:)) = 0.5.*(var(ic,EL2NOD(3,:)) + var(ic,EL2NOD(4,:)));
    var(ic,EL2NOD( 8,:)) = 0.5.*(var(ic,EL2NOD(1,:)) + var(ic,EL2NOD(4,:)));
    var(ic,EL2NOD( 9,:)) = 0.5.*(var(ic,EL2NOD(1,:)) + var(ic,EL2NOD(3,:)));
    var(ic,EL2NOD(10,:)) = 0.5.*(var(ic,EL2NOD(2,:)) + var(ic,EL2NOD(4,:)));
end

end % END OF SUBFUNCTION interp_edge_node_vals

% #########################################################################