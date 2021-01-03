function SAVE_GPLATES_VEL_BCS_TO_MAT_FILE(OPTS)
% Usage: SAVE_GPLATES_VEL_BCS_TO_MAT_FILE(OPTS)
%
% Purpose: 
%   Save GPlates velocity data to a mat file.
%
% Input:
%   OPTS : [structure] : Structure containing info for input and output
%                        files 
%
% Output:
%
% JMT Nov 2017

pdir = pwd;
cd('..');
pwd_main = pwd;
% addpath([pwd_main '/mfiles_M3TET']);
addpath([pwd_main '/mfiles_SPH']);
cd(pdir);

if nargin == 0
    
    OPTS.input_dir = ['/Users/jorge/Tests/SPH_MESH/MESHES_FOR_SOUTH_ATLANTIC_TESTS/' ...
                      'Rectangular_refined_region/n37k_l_ref_100km_l_coarse_2000km/' ... 
                      'n37k_nmg3_Vel_BCs_130_100_Myr'];        % folder where the GPlates data is
    OPTS.outdir    = [pwd_main '/SETUP_TEST/GPLATES/Vel_BCs']; % folder to save the data
    OPTS.name      = 'Vel_BCs_n37k_nmg3';                      % file name to save the data
end

if ismac || isunix
    OPTS.input_dir = change_slash_for_mac(OPTS.input_dir);
    OPTS.outdir    = change_slash_for_mac(OPTS.outdir);
elseif ispc
    OPTS.input_dir = change_slash_for_windows(OPTS.input_dir);
    OPTS.outdir    = change_slash_for_windows(OPTS.outdir);
end

% READ GPLATES VELOCITY MODEL
% ===========================
cd(OPTS.input_dir);
%--------------------------------------------------------------------------
% Check V_BCs file format and number of files
%--------------------------------------------------------------------------
Vel_BCs_info = dir();                % get info about the files in this directory
sample_file  = Vel_BCs_info(4).name; % select 1st Vel_BCs file to check the extension
if strcmp(sample_file(end-2:end),'.xy')
    system('./xy2txt.sh'); % change .xy extension to .txt extension
elseif strcmp(sample_file(end-3:end),'.txt')
    % V_BCs files are already in .txt format
else
    error('wrong extension for V_BCs files. It should be .xy or .txt')
end
Vel_BCs_info = dir('*.txt');
num_files    = numel(Vel_BCs_info); % number of Vel_BCs files

%--------------------------------------------------------------------------
% Create structure for V_BCs and load the files
%--------------------------------------------------------------------------
V_BCs(1:num_files) = struct('time'         ,0, ...
                            'Uth_Uph'      ,0, ...
                            'latlon'       ,0, ...
                            'PlateID'      ,' ', ...
                            'plateid2name' ,' ');
t      = str2double(Vel_BCs_info(num_files).name(1:3)); % initial time (Myr)
t_step = 1; % time step for Vel_BCs files (Myr)
time_VBCs_files = zeros(size(V_BCs,2),1); % vector with the time for each V_BCs file
fprintf('\n Loading GPlates velocity BCs...');
for i = 1:num_files
    time_VBCs_files(i) = t;
    V_BCs(i).time              = t;
    [V_BCs(i).Uth_Uph,V_BCs(i).latlon,V_BCs(i).PlateID,V_BCs(i).plateid2name] = ...
        readVelGPlates([num2str(t),'Ma.txt']);
    t = t - t_step;
end
fprintf(' done.\n\n');
cd(pdir);

clear Vel_BCs_info i numfiles sample_file t t_step VelGPlatesFileName

%--------------------------------------------------------------------------
% Save data
%--------------------------------------------------------------------------
save([OPTS.outdir OPTS.name],'V_BCs','time_VBCs_files');

end % END OF FUNCTION SAVE_GPLATES_VEL_BCS_TO_MAT_FILE

%##########################################################################
%                              SUB-FUNCTIONS
%##########################################################################

function path2data = change_slash_for_windows(path2data)

% '\' can cause problems on unix system; however matlab can handle '/' on
% all architectures
path2data(path2data=='/') = '\';
% make sure there's a slash at the end of the path
if path2data(end)~='\'
    path2data = [path2data '\'];
end

end % END OF SUBFUNCTION change_slash_for_windows

%##########################################################################

function path2data = change_slash_for_mac(path2data)

% '\' can cause problems on unix system; however matlab can handle '/' on
% all architectures
path2data(path2data=='\') = '/';
% make sure there's a slash at the end of the path
if path2data(end)~='/'
    path2data = [path2data '/'];
end

end % END OF SUBFUNCTION change_slash_for_mac