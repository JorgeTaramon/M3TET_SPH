function CREATE_SUBDOMAINS_REF_REGION(OPTS)
% Usage: CREATE_SUBDOMAINS_REF_REGION(OPTS)
%
% Purpose: 
%   Create subdomains of the output data. By default it creates 3
%   subdomains (refined, transition and coarse region) by setting 
%   nsd_ref = 1, but can create more subdomains in the refined region.
%   The function also creates .szplt files to optimize the loading time in
%   Tecplot.
%
% Input:
%   OPTS : [structure] : Structure containing info for creating the
%                        subdomains
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
% (2) specify total number of subdomains in the refined region
%     empty ([]) will find out looking in outdir
nsd_ref = 1; %[];
% (3) specify output directory (folder where data is saved)
outdir2         = '/Users/jorge/Tests/SPH/Trash_01/_Subdomains';
SETTINGS.outdir = outdir2;
create_output_folder(SETTINGS); clear SETTINGS

% Overwrite defaults if different values are specified in structure OPTS
if nargin>0
    if isfield(OPTS,'outdir')
        outdir = OPTS.outdir;
    end
    if isfield(OPTS,'nsd')
        nsd_ref = OPTS.nsd;
    end
end

if numlabs==1
    FigNo_split   = 0;
else
    FigNo_split   = 0;
end

%==========================================================================
% LOAD MESH DATA
%==========================================================================
outdir  = [outdir filesep];
fprintf(' Loading mesh data...');
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
% LOAD VARIABLES
%==========================================================================
GCOORD          = MESH.GCOORD;
GCOORD_SPH      = MESH.GCOORD_SPH;
EL2NOD          = MESH.EL2NOD{1};
nel             = MESH.nel;
nnodel          = size(EL2NOD,1);
% ELEMENTS INSIDE REFINED REGION
els_ref         = MESH.els_ref; 
% ELEMENTS INSIDE TRANSITION REGION
deg2rad         = pi/180;
r_surf          = MESH.r_surf;
theta0          = MESH.theta0; % colatitude (degrees) of the point around which the refined and transition zones are defined
phi0            = MESH.phi0;   % longitude (degrees) of the point around which the refined and transition zones are defined
d_tran          = MESH.d_tran;  % transition zone depth (km), upper mantle - lowe mantle transition
w_tran_deg      = MESH.w_tran/(deg2rad*r_surf); % width of transition zone in degrees (North-South)
theta_tran_n    = theta0 - w_tran_deg/2; % colatitude of the northern boundary in the transition zone
theta_tran_s    = theta0 + w_tran_deg/2; % colatitude of the southern boundary in the transition zone
l_tran_deg      = MESH.l_tran/(deg2rad*r_surf); % length of transition zone in degrees (East-West)
phi_tran_e      = phi0   + l_tran_deg/2; % longitude of the eastern boundary in the transition zone
phi_tran_w      = phi0   - l_tran_deg/2; % longitude of the western boundary in the transition zone
% Compute the barycenter of each element (it determines the element position) in order to know which element are inside the transition region
x_bary          = (GCOORD(1,EL2NOD(1,:)) + GCOORD(1,EL2NOD(2,:)) + GCOORD(1,EL2NOD(3,:)) + GCOORD(1,EL2NOD(4,:)))/4;
y_bary          = (GCOORD(2,EL2NOD(1,:)) + GCOORD(2,EL2NOD(2,:)) + GCOORD(2,EL2NOD(3,:)) + GCOORD(2,EL2NOD(4,:)))/4;
z_bary          = (GCOORD(3,EL2NOD(1,:)) + GCOORD(3,EL2NOD(2,:)) + GCOORD(3,EL2NOD(3,:)) + GCOORD(3,EL2NOD(4,:)))/4;
GCOORD_bary     = [x_bary; y_bary; z_bary];
GCOORD_bary_SPH = cartesian2spherical(GCOORD_bary);
theta           = GCOORD_bary_SPH(1,:)/deg2rad; % theta in degrees
phi             = GCOORD_bary_SPH(2,:)/deg2rad; % phi in degrees
r               = GCOORD_bary_SPH(3,:);
els_tran        = 1:size(EL2NOD,2);
% Take those elements that are inside the transition zone
els_tran        = els_tran(theta >= theta_tran_n  & theta <= theta_tran_s & ...
                           phi   >= phi_tran_w    & phi   <= phi_tran_e & ...
                           r     >= r_surf-d_tran);
% Remove those elements that are inside the refined zone
els_tran        = els_tran(~ismember(els_tran,els_ref));
% ELEMENTS INSIDE COARSE REGION (rest of the mesh)
els_coarse      = 1:size(EL2NOD,2);
% Remove those elements that are inside the refined zone
els_coarse      = els_coarse(~ismember(els_coarse,union(els_ref,els_tran)));
cd(outdir);
Data_info = dir('Lab01x01_Data_*');
num_files = numel(Data_info); % number of Data_info files
fprintf(' done\n');
clear MESH
%==========================================================================
% SPLIT DATA FOR EACH SUBDOMAIN
%==========================================================================
if nsd_ref == 1
    % Split the mesh in 3 subdomains: 
    %   sd1 for the refined region
    %   sd2 for the transition region
    %   sd3 for the coarse region 
    nsd         = nsd_ref + 2;
    el2sd       = zeros(1,size(EL2NOD,2));
    el2sd(ismember((1:nel),els_ref))    = 1*ones(1,size(els_ref,2));
    el2sd(ismember((1:nel),els_tran))   = 2*ones(1,size(els_tran,2));
    el2sd(ismember((1:nel),els_coarse)) = 3*ones(1,size(els_coarse,2));
    
    % Split data for each subdomain
    for isd = 1:nsd
        EL2NOD_SDi         = EL2NOD(:,el2sd==isd);  % rows of "nodes" with elements of SD "isd"
        [nod_SDi2D,~,ind]  = unique(EL2NOD_SDi(:)); % unique list of the nodes in SD "isd"
        nel_SD             = size(EL2NOD_SDi,2);    % number of elements in SD "isd"
        EL2NOD_SDi_P       = EL2NOD(1:4,el2sd==isd);
        [nod_SDi2D_P,~,~]  = unique(EL2NOD_SDi_P(:)); % unique list of the nodes in SD "isd" for P connectivity
        
        % Store subdomain MESH
        MESH.GCOORD = GCOORD(:,nod_SDi2D);                  % coordinates of all nodes in SD "isd
        MESH.EL2NOD = {reshape(uint32(ind),nnodel,nel_SD)}; % element connectivity in SD "isd"
        MESH.nnod   = size(MESH.GCOORD,2);
        prefix      = ['Lab' num2str_d(nsd,2) 'x' num2str_d(isd,2)];
        save([outdir2 '/' prefix '_Mesh'],'MESH');
        
        % Store the subdomain SETTINGS
        load('Lab01x01_SETTINGS.mat')
        save([outdir2 '/' prefix '_SETTINGS'],'SETTINGS','PHYSICS','NUMSCALE');
        
        % Load data for each time step
        fprintf('\n Saving subdomain %d data for each time step (might take a while)...',isd);
        tic
        for i = 1:num_files
            % Store subdomain VAR and time
            load(Data_info(i).name)
            VAR.Ux      = VAR.Ux(nod_SDi2D);
            VAR.Uy      = VAR.Uy(nod_SDi2D);
            VAR.Uz      = VAR.Uz(nod_SDi2D);
            VAR.Uth     = VAR.Uth(nod_SDi2D);
            VAR.Uph     = VAR.Uph(nod_SDi2D);
            VAR.Ur      = VAR.Ur(nod_SDi2D);
            VAR.P       = VAR.P(nod_SDi2D_P);
            VAR.T       = VAR.T(nod_SDi2D);
            VAR.Dens    = VAR.Dens(nod_SDi2D);
            VAR.Visc    = VAR.Visc(nod_SDi2D);
            VAR.Plate   = VAR.Plate(nod_SDi2D);
%             VAR.sigma_n = VAR.sigma_n(nod_SDi2D);
%             VAR.h_dyn   = VAR.h_dyn(nod_SDi2D);
            VAR.h_iso   = VAR.h_iso(nod_SDi2D);
            VAR.h_plume = VAR.h_plume(nod_SDi2D);
%             VAR.M       = VAR.M(nod_SDi2D);
%             VAR.Dpl     = VAR.Dpl(nod_SDi2D);
%             VAR.Vol     = VAR.Vol(nod_SDi2D);
            suffix      = ['_' num2str_d(i-1,4)];
            save([outdir2 '/' prefix '_Data' suffix],'VAR','time');
        end
        fprintf('done.\n');
        toc
    end
    cd(pdir)
    % Plot subdomains
    if FigNo_split
        fprintf('\n Drawing subdomain configuration (might take a few sec)...');
        plot_subdomains(GCOORD,EL2NOD,el2sd,nsd,FigNo_split); % *SUBFUNCTION*
        fprintf('done.\n');
    end
    
elseif nsd_ref > 1 
    % Split the mesh in nsd_ref + 2 subdomains: 
    %   1 to nsd_ref for the refined region
    %   nsd_ref + 1 for the transition region
    %   nsd_ref + 2 for the coarse region 
    cd(pdir)
    EL2NOD_D   = EL2NOD(:,els_ref);
    COMM.nsd   = nsd_ref;
    % Define a percentage of shared nodes that you are willing to "pay" for one less communication cycle. Use values between 5 and ~20.
    % INFO: - A smaller value will try to reduce the number of shared nodes on subdomain boundaries, but more connected 
    %         subdomains will be generated (more communication cycles).
    %       - A higher value does the opposite: It tries to minimize the communication cycles but there qill be more shared nodes.
    fewer_NB_wght = 5;
    el2sd_ref     = el2sd_bisect_3d_sph(GCOORD,GCOORD_SPH,EL2NOD_D,COMM,fewer_NB_wght,1);
    
    el2sd                               = zeros(1,size(EL2NOD,2));
    el2sd(ismember((1:nel),els_ref))    = el2sd_ref';
    el2sd(ismember((1:nel),els_tran))   = (double(max(el2sd_ref)) + 1)*ones(1,size(els_tran,2));
    el2sd(ismember((1:nel),els_coarse)) = (double(max(el2sd_ref)) + 2)*ones(1,size(els_coarse,2));
    
    nsd = double(max(el2sd));
    cd(outdir);
    % Split data for each subdomain
    for isd = 1:nsd
        EL2NOD_SDi         = EL2NOD(:,el2sd==isd);  % rows of "nodes" with elements of SD "isd"
        [nod_SDi2D,~,ind]  = unique(EL2NOD_SDi(:)); % unique list of the nodes in SD "isd"
        nel_SD             = size(EL2NOD_SDi,2);    % number of elements in SD "isd"
        EL2NOD_SDi_P       = EL2NOD(1:4,el2sd==isd);
        [nod_SDi2D_P,~,~]  = unique(EL2NOD_SDi_P(:)); % unique list of the nodes in SD "isd" for P connectivity
        
        % Store subdomain MESH
        MESH.GCOORD = GCOORD(:,nod_SDi2D);                  % coordinates of all nodes in SD "isd
        MESH.EL2NOD = {reshape(uint32(ind),nnodel,nel_SD)}; % element connectivity in SD "isd"
        MESH.nnod   = size(MESH.GCOORD,2);
        prefix      = ['Lab' num2str_d(nsd,2) 'x' num2str_d(isd,2)];
        save([outdir2 '/' prefix '_Mesh'],'MESH');
        
        % Store the subdomain SETTINGS
        load('Lab01x01_SETTINGS.mat')
        save([outdir2 '/' prefix '_SETTINGS'],'SETTINGS','PHYSICS','NUMSCALE');
        
        % Load data for each time step
        fprintf('\n Saving subdomain %d data for each time step (might take a while)...',isd);
        tic
        for i = 1:num_files
            % Store subdomain VAR and time
            load(Data_info(i).name)
            VAR.Ux      = VAR.Ux(nod_SDi2D);
            VAR.Uy      = VAR.Uy(nod_SDi2D);
            VAR.Uz      = VAR.Uz(nod_SDi2D);
            VAR.Uth     = VAR.Uth(nod_SDi2D);
            VAR.Uph     = VAR.Uph(nod_SDi2D);
            VAR.Ur      = VAR.Ur(nod_SDi2D);
            VAR.P       = VAR.P(nod_SDi2D_P);
            VAR.T       = VAR.T(nod_SDi2D);
            VAR.Dens    = VAR.Dens(nod_SDi2D);
            VAR.Visc    = VAR.Visc(nod_SDi2D);
            VAR.Plate   = VAR.Plate(nod_SDi2D);
%             VAR.sigma_n = VAR.sigma_n(nod_SDi2D);
%             VAR.h_dyn   = VAR.h_dyn(nod_SDi2D);
            VAR.h_iso   = VAR.h_iso(nod_SDi2D);
            VAR.h_plume = VAR.h_plume(nod_SDi2D);
            VAR.M       = VAR.M(nod_SDi2D);
            VAR.Dpl     = VAR.Dpl(nod_SDi2D);
            VAR.Vol     = VAR.Vol(nod_SDi2D);
            suffix      = ['_' num2str_d(i-1,4)];
            save([outdir2 '/' prefix '_Data' suffix],'VAR','time');
        end
        fprintf('done.\n');
        toc
    end
    cd(pdir)
    % Plot subdomains
    if FigNo_split
        fprintf('\n Drawing subdomain configuration (might take a few sec)...');
        plot_subdomains(GCOORD,EL2NOD,el2sd,nsd,FigNo_split); % *SUBFUNCTION*
        fprintf('done.\n');
    end
end

%==========================================================================
% CREATE TECPLOT DATA (.dat)
%==========================================================================
% (1) specify output directory (folder where data is located)
OPTS.outdir = outdir2;
% (2) specify which output files are converted
%     (number OR 'all' OR 'last') , e.g. OPTS.plotfiles = [0 113];
OPTS.plotfiles = 'all';
% (3) specify total number of subdomains
%     empty ([]) will find out looking in outdir
OPTS.nsd = nsd; %[];
% (4) specify which subdomain data is converted
%     empty ([]) will convert all subdomains
OPTS.indSD = [];
% (5) specify which variables are converted
OPTS.var_choice = 'all';
% (6) specify suffix of Tecplot file name
%     Tecplot file will be named ['Lab32x24' suffix]
OPTS.suffix = '_Tec360.dat';
BUILD_TEC360FILE_3D(OPTS)

%==========================================================================
% REMOVE .mat FILES
%==========================================================================
cd(outdir2);
delete *.mat
Data_info = dir('*.dat');
cd(pdir)

%==========================================================================
% CONVERT .DAT FILES TO .SZPLT FILES TO REDUCE THE LOADING TIME IN TECPLOT
%==========================================================================
num_files = numel(Data_info); % number of Data_info files
for i=1:num_files
    fid = fopen('convert2szplt.sh','w');
    fprintf(fid,'#!/bin/sh \n');
    fprintf(fid,'TEC=/Applications \n');
    fprintf(fid,'tec360(){ \n');
    fprintf(fid,'$TEC/"Tecplot 360 EX 2016 R2"/"Tecplot 360 EX 2016 R2.app"/Contents/MacOS/"Tecplot 360 EX 2016 R2" "$@" \n');
    fprintf(fid,'} \n');
    fprintf(fid,'export -f tec360 \n');
    fprintf(fid,'tec360 -b -p convert2szplt.mcr %s',outdir2);
    fprintf(fid,'/%s \n',Data_info(i).name);
    fclose(fid);
    system('bash ./convert2szplt.sh')
end

fid = fopen('dat_szplt2szplt.sh','w');
fprintf(fid,'#!/bin/sh \n');
fprintf(fid,'cd %s \n',outdir2);
fprintf(fid,'for old in *.dat.szplt \n');
fprintf(fid,'do mv $old `basename $old .dat.szplt`.szplt \n');
fprintf(fid,'done \n');
fclose(fid);
system('bash ./dat_szplt2szplt.sh'); % change .dat.szplt extension to .szplt extension
cd(outdir2);
delete *.dat
cd(pdir)
delete batch.log
delete convert2szplt.sh
delete dat_szplt2szplt.sh
end

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function plot_subdomains(GCOORD,EL2NOD,el2sd,nsd,FigNo_split)

figure(FigNo_split);clf
cmap = jet(nsd);

GCOORD(1,:) = GCOORD(1,:) - 5000;
GCOORD(2,:) = GCOORD(2,:) - 5000;
for isd=1:nsd
    simpplot(GCOORD',EL2NOD(:,el2sd==isd)','p(:,1)<-5000 & p(:,2)<-5000',cmap(isd,:),cmap(isd,:));
    hold on
end
GCOORD(2,:) = GCOORD(2,:) + 10000;
for isd=1:nsd
    simpplot(GCOORD',EL2NOD(:,el2sd==isd)','p(:,1)<-5000 & p(:,2)>5000',cmap(isd,:),cmap(isd,:));
    hold on
end
GCOORD(1,:) = GCOORD(1,:) + 10000;
for isd=1:nsd
    simpplot(GCOORD',EL2NOD(:,el2sd==isd)','p(:,1)>5000 & p(:,2)>5000',cmap(isd,:),cmap(isd,:));
    hold on
end
GCOORD(2,:) = GCOORD(2,:) - 10000;
for isd=1:nsd
    simpplot(GCOORD',EL2NOD(:,el2sd==isd)','p(:,1)>5000 & p(:,2)<-5000',cmap(isd,:),cmap(isd,:));
    hold on
end
view(150,40); xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');

end % END OF SUBFUNCTION plot_subdomains