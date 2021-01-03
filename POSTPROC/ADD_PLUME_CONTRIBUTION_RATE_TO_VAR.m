function ADD_PLUME_CONTRIBUTION_RATE_TO_VAR(OPTS)
% Usage: ADD_PLUME_CONTRIBUTION_RATE_TO_VAR(OPTS)
%
% Purpose: 
%   Add plume contribution rate (uplift rate) to output data:
%   
%                        h_plume(i+1) - h_plume (i)
%        h_dot_right = ------------------------------
%                              t(i+1) - t(i)
%
%                        h_plume(i) - h_plume (i-1)
%        h_dot_left = ------------------------------
%                              t(i) - t(i-1)
%   
%   where i is the time step. Then,
%   
%                 h_dot_right + h_dot_left
%        h_dot = ------------------------------
%                             2
%
%   IMPORTANT: This routine can only be applied after using
%   ADD_PLUME_CONTRIBUTION_TO_VAR, where h_plume is computed
%
% Input:
%   OPTS : [structure] : Structure containing info for output file
%
% Output:
%
% JMT Oct 2018

pdir = pwd;
cd('..');
addpath([pwd '/mfiles_M3TET']);
addpath([pwd '/mfiles_SPH']);
cd(pdir);


% Specify output directory for the plume model(folder where data is located)
outdir_plume       = '/Users/jorge/Tests/SPH/TESTS_SOUTH_ATLANTIC_GPLATES/L1F20';
outdir_plume       = [outdir_plume filesep];
outdir_plume_extra = '/Users/jorge/Tests/SPH/TESTS_SOUTH_ATLANTIC_GPLATES/L1F20/extra_files';
outdir_plume_extra = [outdir_plume_extra filesep];

% Overwrite defaults if different values are specified in structure OPTS
if nargin>0
    if isfield(OPTS,'outdir')
        outdir_plume = OPTS.outdir;
    end
end
% Specify output directory where data will be saved
SETTINGS.outdir = '/Users/jorge/Tests/SPH/L1F20_uplift_rate/';
create_output_folder(SETTINGS);
% Specify the name of the variables for .plt files
varnames = {'Ux' 'Uy' 'Uz' 'Uth' 'Uph' 'Ur' 'P' 'T' 'Dens' 'Visc' 'Plate' 'h_iso' 'h_plume' 'h_dot'};

%==========================================================================
% LOAD MESH DATA
%==========================================================================
fprintf(' Loading mesh data and settings...');
try
    data  = load([outdir_plume 'Lab01x01_Mesh.mat']);
catch
    error(' Could no open file "%s"\n',[outdir_plume 'MESH.mat']);
end
if isfield(data,'Mesh')
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
    data2  = load([outdir_plume 'Lab01x01_SETTINGS.mat']);
catch
    error(' Could no open file "%s"\n',[outdir_plume 'SETTINGS.mat']);
end
NUMSCALE = data2.NUMSCALE;
clear data2
fprintf(' done\n');

%==========================================================================
% COMPUTE PLUME CONTRIBUTION RATE (UPLIFT RATE) (0 Myr - 20 Myr)
%==========================================================================
cd(outdir_plume);
Data_info_plume       = dir('Lab01x01_Data_*');
cd(outdir_plume_extra);
Data_info_plume_extra = dir('Lab01x01_Data_*');
cd(outdir_plume);

time_Myr   = [0  2  4  6  8  10  12  14  16  18  20]; % (Myr)
Data_file  = [1 21 42 62 82 102 122 142 162 182 202];
num_file   = Data_file + 1; % we have to add 1 because the first file is 0.

for i = 1:size(num_file,2)
   if num_file(i) == num_file(1)
       load(Data_info_plume(num_file(i)+1).name) % load VAR in file i+1
       h_plume_i_plus_1 = VAR.h_plume;
       load(Data_info_plume(num_file(i)).name)   % load VAR in file i
       % Compute uplift rate (m/Myr)
       VAR.h_dot = 1000*(h_plume_i_plus_1 - VAR.h_plume)/0.1;
       suffix    = ['_' num2str_d(Data_file(i),4)];
       save([SETTINGS.outdir 'Lab01x01_Data' suffix],'dt','istep','time','VAR');
       write_tecplot_data_3d(MESH,VAR,SETTINGS,NUMSCALE,time,Data_file(i),Data_file(i),varnames);
       fprintf(1,'================================================================\n');
       fprintf('Time: %2i Myr', time_Myr(i))
       fprintf(1,'\n')
       fprintf('Max h_dot: %5.2f m/Myr', max(VAR.h_dot))
       fprintf(1,'\n')
       fprintf('Min h_dot: %5.2f m/Myr', min(VAR.h_dot))
       fprintf(1,'\n')
   elseif num_file(i) == num_file(end)
       cd(outdir_plume_extra);
       load(Data_info_plume_extra(1).name) % load VAR in file 1 in outdir_plume_extra
       h_plume_i_plus_1 = VAR.h_plume;
       cd(outdir_plume);
       load(Data_info_plume(num_file(i)-1).name) % load VAR in file i-1
       h_plume_i_minus_1 = VAR.h_plume;
       load(Data_info_plume(num_file(i)).name)   % load VAR in file i
       % Compute uplift rate (m/Myr)
       VAR.h_dot = 1000*((h_plume_i_plus_1 - VAR.h_plume)/0.1 + (VAR.h_plume - h_plume_i_minus_1)/0.1)/2;
       suffix    = ['_' num2str_d(Data_file(i),4)];
       save([SETTINGS.outdir 'Lab01x01_Data' suffix],'dt','istep','time','VAR');
       write_tecplot_data_3d(MESH,VAR,SETTINGS,NUMSCALE,time,Data_file(i),Data_file(i),varnames);
       fprintf(1,'================================================================\n');
       fprintf('Time: %2i Myr', time_Myr(i))
       fprintf(1,'\n')
       fprintf('Max h_dot: %5.2f m/Myr', max(VAR.h_dot))
       fprintf(1,'\n')
       fprintf('Min h_dot: %5.2f m/Myr', min(VAR.h_dot))
       fprintf(1,'\n')
   else
       load(Data_info_plume(num_file(i)+1).name) % load VAR in file i+1
       h_plume_i_plus_1 = VAR.h_plume;
       load(Data_info_plume(num_file(i)-1).name) % load VAR in file i-1
       h_plume_i_minus_1 = VAR.h_plume;
       load(Data_info_plume(num_file(i)).name)   % load VAR in file i
       % Compute uplift rate (m/Myr)
       VAR.h_dot = 1000*((h_plume_i_plus_1 - VAR.h_plume)/0.1 + (VAR.h_plume - h_plume_i_minus_1)/0.1)/2;
       suffix    = ['_' num2str_d(Data_file(i),4)];
       save([SETTINGS.outdir 'Lab01x01_Data' suffix],'dt','istep','time','VAR');
       write_tecplot_data_3d(MESH,VAR,SETTINGS,NUMSCALE,time,Data_file(i),Data_file(i),varnames);
       fprintf(1,'================================================================\n');
       fprintf('Time: %2i Myr', time_Myr(i))
       fprintf(1,'\n')
       fprintf('Max h_dot: %5.2f m/Myr', max(VAR.h_dot))
       fprintf(1,'\n')
       fprintf('Min h_dot: %5.2f m/Myr', min(VAR.h_dot))
       fprintf(1,'\n')
   end
end

%================================================================================
% COMPUTE PLUME CONTRIBUTION RATE (UPLIFT RATE) FOR EXTRA FILES (20 Myr - 30 Myr)
%================================================================================
cd(outdir_plume_extra);
Data_info_plume_extra = dir('Lab01x01_Data_*');

time_extra_Myr  = [ 22  24  26  28  30]; % (Myr)
Data_file_extra = [222 242 262 282 302];
num_file_extra  = [ 20  40  60  80 100];

for i = 1:size(Data_file_extra,2)
   if Data_file_extra(i) == Data_file_extra(end)
       load(Data_info_plume_extra(num_file_extra(i)-1).name) % load VAR in file i-1
       h_plume_i_minus_1 = VAR.h_plume;
       load(Data_info_plume_extra(num_file_extra(i)).name)   % load VAR in file i
       % Compute uplift rate (m/Myr)
       VAR.h_dot = 1000*(VAR.h_plume - h_plume_i_minus_1)/0.1;
       suffix    = ['_' num2str_d(Data_file_extra(i),4)];
       save([SETTINGS.outdir 'Lab01x01_Data' suffix],'dt','istep','time','VAR');
       write_tecplot_data_3d(MESH,VAR,SETTINGS,NUMSCALE,time,Data_file_extra(i),Data_file_extra(i),varnames);
       fprintf(1,'================================================================\n');
       fprintf('Time: %2i Myr', time_extra_Myr(i))
       fprintf(1,'\n')
       fprintf('Max h_dot: %5.2f m/Myr', max(VAR.h_dot))
       fprintf(1,'\n')
       fprintf('Min h_dot: %5.2f m/Myr', min(VAR.h_dot))
       fprintf(1,'\n')
   else
       load(Data_info_plume_extra(num_file_extra(i)+1).name) % load VAR in file i+1
       h_plume_i_plus_1 = VAR.h_plume;
       load(Data_info_plume_extra(num_file_extra(i)-1).name) % load VAR in file i-1
       h_plume_i_minus_1 = VAR.h_plume;
       load(Data_info_plume_extra(num_file_extra(i)).name)   % load VAR in file i
       % Compute uplift rate (m/Myr)
       VAR.h_dot = 1000*((h_plume_i_plus_1 - VAR.h_plume)/0.1 + (VAR.h_plume - h_plume_i_minus_1)/0.1)/2;
       suffix    = ['_' num2str_d(Data_file_extra(i),4)];
       save([SETTINGS.outdir 'Lab01x01_Data' suffix],'dt','istep','time','VAR');
       write_tecplot_data_3d(MESH,VAR,SETTINGS,NUMSCALE,time,Data_file_extra(i),Data_file_extra(i),varnames);
       fprintf(1,'================================================================\n');
       fprintf('Time: %2i Myr', time_extra_Myr(i))
       fprintf(1,'\n')
       fprintf('Max h_dot: %5.2f m/Myr', max(VAR.h_dot))
       fprintf(1,'\n')
       fprintf('Min h_dot: %5.2f m/Myr', min(VAR.h_dot))
       fprintf(1,'\n')
   end  
end
end % END OF FUNCTION ADD_PLUME_CONTRIBUTION_RATE_TO_VAR