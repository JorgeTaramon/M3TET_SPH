function CREATE_CSV_FILES_FOR_COASTLINES(OPTS)
% Usage: CREATE_CSV_FILES_FOR_COASTLINES(OPTS)
%
% Purpose:
%   Create csv file with the geometry of the present day coastlines
%   reconstructed every 2 Ma from 130 to 100 Ma in order to be read by
%   Tecplot. Steps:
%   - Load mesh (where rotation matrices are)
%   - Read the coastline files (lat,lon).
%   - Transform them into Cartesian coordinates and rotate them to place
%     near the equator.
%   - Create csv files so that coastlines can be loaded in Tecplot.
%   - Repeat for every 2Ma
%
% Input:
%   OPTS : [structure] : Structure containing info for output file
%
% Output:
%
% JMT Oct 2017

cd('..');
addpath([pwd '/mfiles_M3TET']);
addpath([pwd '/mfiles_SPH']);

% specify output directory (folder where data is going to be saved)
% outdir = '/Users/jorge/Tests/SPH/TESTS_SOUTH_ATLANTIC_GPLATES/1.Initial_tests/n64k_nmg3(=4M)_rectangular_ref_region_130_100_Ma/Plume(48.3d)';
outdir = '/Users/jorge/Tests/SPH/TESTS_SOUTH_ATLANTIC_GPLATES';
% Overwrite defaults if different values are specified in structure OPTS
if nargin>0
    if isfield(OPTS,'outdir')
        outdir = OPTS.outdir;
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


cd SETUP_TEST/GPLATES/COASTLINES; % directory where the coastlines files are
pdir = pwd;

D = dir; % info about the current directory
for i = 4:length(D)
    current_folder = D(i).name; % Get the current subdirectory name
    
    cd(current_folder)
    
    %=================================================================================
    % CHECK COASTLINES FILE FORMAT AND NUMBER OF FILES
    %=================================================================================
    coastlines_info = dir();                   % get info about the files in this directory
    sample_file     = coastlines_info(4).name; % select 1st coastline file to check the extension
    if strcmp(sample_file(end-2:end),'.xy')
        system('./xy2txt.sh'); % change .xy extension to .txt extension
    elseif strcmp(sample_file(end-3:end),'.txt')
        % Coastlines files are already in .txt format
    else
        error('wrong extension for COASTLINES files. It should be .xy or .txt')
    end
    coastlines_info = dir('*.txt');
    num_files    = numel(coastlines_info); % number of coastlines files
    
    %=================================================================================
    % LOAD FILES, TRANSFORM (lat,lon) INTO (x,y,z) AND ROTATE THEM FROM GPLATES FRAME
    %=================================================================================
    for j = 1:num_files
        fid     = fopen(coastlines_info(j).name,'r'); % open input file
        tmp     = fscanf(fid, '%f');
        lat_lon = reshape(tmp,2,[])';
        fclose(fid); % close input file
        
        theta                 = (90 - lat_lon(:,1))*pi/180; % colatitude in radians
        phi                   = lat_lon(:,2)*pi/180;        % longitude in radians
        phi(phi < 0) = phi(phi < 0) + 2*pi;                 % values between 0 and 2pi
        GCOORD_SPH_coastlines = [theta phi repmat(MESH.r_surf,size(lat_lon,1),1)]';

        GCOORD_coastlines     = spherical2cartesian(GCOORD_SPH_coastlines);
        GCOORD_coastlines_rot = (MESH.RR2 * MESH.RR1)' * GCOORD_coastlines; % rotate coordinates from GPlates frame
        GCOORD_coastlines_rot = GCOORD_coastlines_rot';
        
        %=============================================================================
        % CREATE CSV FILES WITH HEADERS
        %=============================================================================
        headers = {'X','Y','Z'}; % headers to save in the csv file (needs the funtion csvwrite_with_headers.m)
        cd(outdir)
        cd('Coastlines')
        cd(current_folder)
        if j == 1
            csvwrite_with_headers('Africa_1.csv',GCOORD_coastlines_rot,headers)
        elseif j == 2
            csvwrite_with_headers('Africa_2.csv',GCOORD_coastlines_rot,headers)
        elseif j == 3
            csvwrite_with_headers('South_America_1.csv',GCOORD_coastlines_rot,headers)
        elseif j == 4
            csvwrite_with_headers('South_America_2.csv',GCOORD_coastlines_rot,headers)
        elseif j == 5
            csvwrite_with_headers('South_America_3.csv',GCOORD_coastlines_rot,headers)
        elseif j == 6
            csvwrite_with_headers('South_America_4.csv',GCOORD_coastlines_rot,headers)
        elseif j == 7
            csvwrite_with_headers('South_America_5.csv',GCOORD_coastlines_rot,headers)
        elseif j == 8
            csvwrite_with_headers('South_America_6.csv',GCOORD_coastlines_rot,headers)
        elseif j == 9
            csvwrite_with_headers('South_America_7.csv',GCOORD_coastlines_rot,headers)
        end
        cd(pdir)
        cd(current_folder)
    end
    cd('..')
end

end % END OF FUNCTION CREATE_CSV_FILES_FOR_COASTLINES