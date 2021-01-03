function ADD_PLUME_CONTRIBUTION_TO_VAR(OPTS)
% Usage: ADD_PLUME_CONTRIBUTION_TO_VAR(OPTS)
%
% Purpose: 
%   Add plume contribution to output data subtracting mantle contribution
%   to isostatic topography from the plume model minus mantle contribution
%   to isostatic topography from the plume free model.
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


% specify output directory for the plume model(folder where data is located)
outdir_plume = '/Users/jorge/Tests/SPH/TESTS_SOUTH_ATLANTIC_GPLATES_MERGED_CRATONS_8_CONTOURS/L1F15';
% Overwrite defaults if different values are specified in structure OPTS
if nargin>0
    if isfield(OPTS,'outdir')
        outdir_plume = OPTS.outdir;
    end
end
% specify output directory where the data for mantle contribution of the plume free model is located
outdir_no_plume ='/Users/jorge/Tests/SPH/TESTS_SOUTH_ATLANTIC_GPLATES_MERGED_CRATONS_8_CONTOURS/NoPlume/h_iso';

%==========================================================================
% COMPUTE PLUME CONTRIBUTION TO ISOSTATIC TOPOGRAPHY
%==========================================================================
outdir_plume    = [outdir_plume filesep];
outdir_no_plume = [outdir_no_plume filesep];
cd(outdir_plume);
Data_info_plume = dir('Lab01x01_Data_*');
num_files = numel(Data_info_plume); % number of Data_info_plume files
cd(outdir_no_plume);
Data_info_no_plume = dir('h_iso_no_plume_*');

for i = 1:num_files
    cd(outdir_no_plume);
    load(Data_info_no_plume(i).name)
    
    cd(outdir_plume);
    load(Data_info_plume(i).name)
    
    VAR.h_plume = VAR.h_iso - h_iso_no_plume;
    
    suffix = ['_' num2str_d(i-1,4)];
    save([outdir_plume 'Lab01x01_Data' suffix],'dt','istep','time','VAR');
end
cd(pdir)

end % END OF FUNCTION ADD_PLUME_CONTRIBUTION_TO_VAR