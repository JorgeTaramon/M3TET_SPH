function EXTRACT_MANTLE_CONTRIBUTION(OPTS)
% Usage: EXTRACT_MANTLE_CONTRIBUTION(OPTS)
%
% Purpose: 
%   Extract mantle contribution from output data of a plume free model.
%
%   IMPORTANT: This routine can only be applied after using
%   ADD_STRESSES_TO_VAR, where h_iso is computed
%
% Input:
%   OPTS : [structure] : Structure containing info for output file
%
% Output:
%
% JMT Sept 2017


pdir = pwd;
cd('..');
addpath([pwd '/mfiles_M3TET']);
addpath([pwd '/mfiles_SPH']);
cd(pdir);

%=======================================================================================
% EXTRACT h_iso (MANTLE CONTRIBUTION TO ISOSTATIC TOPOGRAPHY) FROM THE PLUME FREE MODEL
%=======================================================================================
% specify output directory for the plume free model(folder where data is located)
outdir_no_plume = '/Users/jorge/Tests/SPH/Trash_01';
% Overwrite defaults if different values are specified in structure OPTS
if nargin>0
    if isfield(OPTS,'outdir')
        outdir_no_plume = OPTS.outdir;
    end
end
% specify output directory (folder where data is saved)
outdir2         ='/Users/jorge/Tests/SPH/Trash_01/h_iso';
SETTINGS.outdir = outdir2;
create_output_folder(SETTINGS); clear SETTINGS

%==========================================================================
% EXTRACT h_iso from VAR
%==========================================================================
outdir    = [outdir_no_plume filesep];
cd(outdir);
Data_info = dir('Lab01x01_Data_*');
num_files = numel(Data_info); % number of Data_info files

for i = 1:num_files
    % Store subdomain VAR and time
    load(Data_info(i).name)
    
    h_iso_no_plume = VAR.h_iso;
    
    suffix = ['_' num2str_d(i-1,4)];
    save([outdir2 '/' 'h_iso_no_plume' suffix],'h_iso_no_plume');
end
cd(pdir)

end % END OF FUNCTION EXTRACT_MANTLE_CONTRIBUTION