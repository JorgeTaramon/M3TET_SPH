function [MESH,COMM] = init_mesh_straight_p(SETTINGS,PHYSICS,NUMSCALE,COMM)
% Usage: [MESH,COMM] = init_mesh_straight_p(SETTINGS,PHYSICS,NUMSCALE,COMM)
%
% Purpose: Load mesh from file, create multigrid levels, store all data in
%          structure MESH. Create structure COMM storing communication
%          infrastructure (also needed if code runs in serial).
%          When parallel: create subdomains
%
% Input:
%   SETTINGS : [structure] : model parameters
%   PHYSICS  : [structure] : physical properties
%   COMM     : [structure] : inter-subdomain communication data
%
% Output:
%   MESH     : [structure] : FE mesh parameters
%   COMM     : [structure] : inter-subdomain communication data
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Mar 2013
% JH May 2014: now prepares workspace for tsearch2 ("WS")
% JH Feb 2016: added comments, cleaned up
%

fidl = SETTINGS.fid_log;

%==========================================================================
% DEFINE WHICH CORNERS AND SURFACES ARE CONNECTED TO WHICH SURFACES
% - nodes at corners have PointIDs 101:199
% - nodes on edges have PointIDs 201:299
% - nodes on surfaces have PointIDs 301:399
% IF "SETTINGS.DB_indices" IS NOT BEEN DEFINED, WE USE THE STANDARD CONFIG
% 
% ON A SIMPLE SPHERE THERE"S ONLY 301 (bottom) and 306 (top)
%==========================================================================
if ~isfield(SETTINGS,'DB_indices')
%     DB_indices    = cell(6,1);
%     % Specify which corners and edges are located on which domain face:
%     %                        |<-- corners -->|<--  edges  -->|face
%     DB_indices{1} = [101 102 103 104 201 202 203 204 301];
%     DB_indices{2} = [101 102 105 106 201 205 206 209 302];
%     DB_indices{3} = [102 103 106 107 202 206 207 210 303];
%     DB_indices{4} = [103 104 107 108 203 207 208 211 304];
%     DB_indices{5} = [101 104 105 108 204 205 208 212 305];
%     DB_indices{6} = [105 106 107 108 209 210 211 212 306];
    % A sphere, however, only has a surface at bottom (301) and top (306)
    DB_indices{1} = 301;
    DB_indices{2} = 306;
else
    DB_indices    = SETTINGS.DB_indices;
end


%==========================================================================
% LOAD COARSEST MULTIGRID MESH
%==========================================================================
nnodel   = 10; % M3TET_SPH needs quadratic tetrahedral elements
[GCOORD,EL2NOD,PointID,PhaseID,r_cmb,r_surf,EMBEDDED_MESH] = ...
    load_mesh_sph(SETTINGS.meshfile,nnodel,fidl,DB_indices);

if SETTINGS.nmg>1
    nnod_coarse_limit = 95000;
    nnod              = size(GCOORD,2);
    if nnod>nnod_coarse_limit
        error('The coarse Multigrid-mesh has more nodes than the current limit (%i>%i)',...
            nnod,nnod_coarse_limit);
    end
end


% %==========================================================================
% % CHECK THAT EDGE NODES ARE CENTERED BETWEEN VERTICES (only straight case)
% %==========================================================================
% nodE = 5:10;
% nodV = [1 2; 2 3; 3 4; 4 1; 1 3; 2 4];
% for i=1:length(nodE)
%     x_test = 0.5*sum(reshape(GCOORD(1,EL2NOD(nodV(i,:),:)),2,[]),1);
%     y_test = 0.5*sum(reshape(GCOORD(2,EL2NOD(nodV(i,:),:)),2,[]),1);
%     z_test = 0.5*sum(reshape(GCOORD(3,EL2NOD(nodV(i,:),:)),2,[]),1);
%     diff   = abs(x_test - GCOORD(1,EL2NOD(nodE(i),:))) + ...
%              abs(y_test - GCOORD(2,EL2NOD(nodE(i),:))) + ...
%              abs(z_test - GCOORD(3,EL2NOD(nodE(i),:)));
%     tol    = 1e-12 * (max(x_test)-min(x_test));
%     ind    = diff>tol;
%     if any(ind)
%         ind = find(ind);
%         if length(ind)==size(EL2NOD,2);
%             error(' This mesh has the wrong connectivity! Node #%i must be between node %i and %i.',nodE,nodV);
%         else
%             error(' Some edge nodes #%i are not centered between corner nodes %i and %i.',nodE,nodV);
%         end
%     end
% end


%==========================================================================
% GENERATE MULTIGRID MESHES
%==========================================================================
if ~isfield(SETTINGS,'nmg')
    error('SETTINGS.nmg is not defined.');
else
    nmg = SETTINGS.nmg;
end

if COMM.nsd>1
    %======================================================================
    % CODE RUNS IN PARALLEL MODE
    % Have to generate subdomains, communication arrays, multigrid levels
    %======================================================================
    error('Not yet coded');
    
    % Save mesh data in structure
    % ===========================
    MESH  = struct('name'         ,SETTINGS.meshfile, ...
                   'GCOORD'       ,GCOORD           , ...
                   'EL2NOD'       ,EL2NOD           , ...
                   'GCOORD_SPH'   ,GCOORD_SPH       , ...
                   'nnodel'       ,nnodel           , ...
                   'nvertx'       ,4                , ...
                   'nnod'         ,size(GCOORD,2)   , ...
                   'nel'          ,size(EL2NOD,2)   , ...
                   'PointID'      ,PointID          , ...
                   'PhaseID'      ,PhaseID          , ...
                   'DB_indices'   ,{DB_indices});
    if nmg>1
        error('it needs to tested if it works well for curved edges meshes') 
        [MESH,COMM] = create_subdomain_multigrid_meshes_3d(MESH,COMM,nmg,fidl);
    else
        error('M3TET_SPH running in parallel requires multigrid (SETTINGS.nmg>1)');
    end

    % Prepare workspace for MUTILS' tsearch2
    % ======================================
    MESH.WS_tsearch2 = calc_tsearch2_WS(MESH.GCOORD_D,MESH.EL2NOD_D);

else
    %======================================================================
    % CODE RUNS IN SERIAL MODE
    % Have to generate only multigrid levels
    %======================================================================
    if nmg>1
        fprintf(fidl,' Creating %1i multigrid level (straight-edges)...',nmg);
        [GCOORD,EL2NOD,PointID,PhaseID,Ic2f] = create_multigrid_meshes_3d...
            (GCOORD,EL2NOD,PointID,PhaseID,DB_indices,nmg);
        fprintf(fidl,'done.\n');
    else
        Ic2f   = [];
        EL2NOD = {uint32(EL2NOD)};
    end
    
    % Compute elements inside the refined region
    if ~isempty(EMBEDDED_MESH.l0_ref)
        [els_ref,points_ref] = compute_els_inside_refined_region(GCOORD,EL2NOD{1},EMBEDDED_MESH,r_surf);
    else
        els_ref    = [];
        points_ref = [];
    end
    
    GCOORD_SPH = cartesian2spherical(GCOORD);
    
    % Save mesh in structure
    % ======================
    MESH = struct('name'       ,SETTINGS.meshfile    , ...
                  'GCOORD'     ,GCOORD               , ...
                  'GCOORD_SPH' ,GCOORD_SPH           , ...
                  'EL2NOD'     ,{EL2NOD}             , ...
                  'els_ref'    ,els_ref              , ...
                  'points_ref' ,points_ref           , ...
                  'r_cmb'      ,r_cmb                , ...
                  'r_surf'     ,r_surf               , ...
                  'theta_cone' ,SETTINGS.theta_cone  , ...
                  'nnodel'     ,nnodel               , ...
                  'nvertx'     ,4                    , ...
                  'nnod'       ,size(GCOORD,2)       , ...
                  'nel'        ,size(EL2NOD{1},2)    , ...
                  'nmg'        ,nmg                  , ...
                  'PointID'    ,PointID              , ...
                  'PhaseID'    ,PhaseID              , ...
                  'DB_indices' ,{DB_indices}         , ...
                  'Ic2f'       ,{Ic2f}               , ...
                  'RR1'        ,EMBEDDED_MESH.RR1    , ...
                  'RR2'        ,EMBEDDED_MESH.RR2    , ...
                  'theta0'     ,EMBEDDED_MESH.theta0 , ...
                  'phi0'       ,EMBEDDED_MESH.phi0   , ...
                  'd_ref'      ,EMBEDDED_MESH.d_ref  , ...
                  'w_ref'      ,EMBEDDED_MESH.w_ref  , ...
                  'l_ref'      ,EMBEDDED_MESH.l_ref);
    
    % Prepare workspace for MUTILS' tsearch2
    % ======================================
    MESH.WS_tsearch2 = calc_tsearch2_ws(MESH.GCOORD,MESH.EL2NOD{nmg});
end

%==========================================================================
% CALCULATE MORE DATA ON THE MESH
%==========================================================================
if nnodel==10
    MESH.nVnod = max(max(MESH.EL2NOD{1}(1:4,:)));
end

for img=1:nmg
    MESH.nnod_mg(img) = max(max(MESH.EL2NOD{img}));
    if nnodel==10
        MESH.nVnod_mg(img) = max(max(MESH.EL2NOD{img}(1:4,:)));
    end
end

% Calculate length of each element
len_el       = calc_tetra_length(MESH.GCOORD,MESH.EL2NOD{1});
MESH.len_el  = [min(len_el) max(len_el)]; % only for information
clear len_el

% Calculate volume of each element
vol_el       = calc_tetra_volume(MESH.GCOORD,MESH.EL2NOD{1});
MESH.vol_el  = [min(vol_el) max(vol_el)]; % only for information

% Calculate volume associated with each node
subs         = tetmesh_p2_to_p1([],MESH.EL2NOD{1});
vol_el       = (1/8) * repmat(vol_el(:)',8,1);
vals    	 = repmat(0.25*vol_el(:)',4,1);
MESH.vol_nod = accumarray(subs(:),vals(:));
clear subs vals vol_el

MESH.xmin   = min(GCOORD(1,:));
MESH.xmax   = max(GCOORD(1,:));
MESH.ymin   = min(GCOORD(2,:));
MESH.ymax   = max(GCOORD(2,:));
MESH.zmin   = min(GCOORD(3,:));
MESH.zmax   = max(GCOORD(3,:));
MESH.Lx     = MESH.xmax-MESH.xmin;
MESH.Ly     = MESH.ymax-MESH.ymin;
MESH.Lz     = MESH.zmax-MESH.zmin;
MESH.r_ftc  = @(inod) sqrt(sum(GCOORD(:,inod).^2,1));

% CALCULATE MATRIX THAT TRANSFORMS BETWEEN CARTESIAN AND SPHERICAL COORDS
switch SETTINGS.plate_model
    case 'WJM'% For WJW plate velocity model
        [MESH.RR_xyz2rEN,MESH.xyzR,MESH.xyzE,MESH.xyzN] = make_rot_matrix(MESH.GCOORD);
    case 'GPlates'% For GPlates/Debug plate velocity model 
        [MESH.RR]                                       = make_rotation_matrix(GCOORD);
         % (same as in the 3d mesh generator)
    case 'Debug'% For Debug
        [MESH.RR]                                       = make_rotation_matrix(GCOORD);
         % (same as in the 3d mesh generator)    
end

if numlabs==1 && COMM.nsd>1 % WE ARE IN TEST MODE
    display_mesh_prop(MESH,COMM,SETTINGS,NUMSCALE);
    error(' Stopping code because one worker was only simulating a parallel mesh initialization.');
end

% % Save mesh data in output folder
% save([SETTINGS.outdir '/' COMM.prefix '_MESH'],'MESH');
% save([SETTINGS.outdir '/' COMM.prefix '_COMM'],'COMM');

end % END OF FUNCTION init_mesh_straight_p

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function [els_ref,points_ref] = compute_els_inside_refined_region(GCOORD,EL2NOD,EMBEDDED_MESH,r_surf)

% COMPUTE EACH ELEMENT BARYCENTER (IT DETERMINES THE ELEMENT POSITION) IN ORDER TO KNOW WHICH ELEMENTS ARE INSIDE THE REFINED REGION 
deg2rad         = pi/180;
theta0          = EMBEDDED_MESH.theta0; % colatitude (degrees) of the point around which the refined and transition zones are defined
phi0            = EMBEDDED_MESH.phi0;   % longitude (degrees) of the point around which the refined and transition zones are defined
d_ref           = EMBEDDED_MESH.d_ref;  % refined zone depth (km)
w_ref_deg       = EMBEDDED_MESH.w_ref/(deg2rad*r_surf);  % width of refined zone in degrees (North-South)
theta_ref_n     = theta0 - w_ref_deg/2; % colatitude of the northern boundary in the refined zone
theta_ref_s     = theta0 + w_ref_deg/2; % colatitude of the southern boundary in the refined zone
l_ref_deg       = EMBEDDED_MESH.l_ref/(deg2rad*r_surf);  % length of refined zone in degrees (East-West)
phi_ref_e       = phi0   + l_ref_deg/2; % longitude of the eastern boundary in the refined zone
phi_ref_w       = phi0   - l_ref_deg/2; % longitude of the western boundary in the refined zone
l0_ref          = EMBEDDED_MESH.l0_ref; % element length inside refined region

x_bary          = (GCOORD(1,EL2NOD(1,:)) + GCOORD(1,EL2NOD(2,:)) + GCOORD(1,EL2NOD(3,:)) + GCOORD(1,EL2NOD(4,:)))/4;
y_bary          = (GCOORD(2,EL2NOD(1,:)) + GCOORD(2,EL2NOD(2,:)) + GCOORD(2,EL2NOD(3,:)) + GCOORD(2,EL2NOD(4,:)))/4;
z_bary          = (GCOORD(3,EL2NOD(1,:)) + GCOORD(3,EL2NOD(2,:)) + GCOORD(3,EL2NOD(3,:)) + GCOORD(3,EL2NOD(4,:)))/4;
GCOORD_bary     = [x_bary; y_bary; z_bary];
GCOORD_bary_SPH = cartesian2spherical(GCOORD_bary);
theta           = GCOORD_bary_SPH(1,:)/deg2rad; % theta in degrees
phi             = GCOORD_bary_SPH(2,:)/deg2rad; % phi in degrees
r               = GCOORD_bary_SPH(3,:);
els_ref         = 1:size(EL2NOD,2);
% Take those elements that are inside the refined zone
els_ref         = els_ref(theta >= theta_ref_n  & theta <= theta_ref_s & ...
                          phi   >=   phi_ref_w  & phi   <=   phi_ref_e & ...
                          r     >= r_surf-d_ref);
% Compute points of the mesh inside the refined region
GCOORD_SPH = cartesian2spherical(GCOORD);
theta      = GCOORD_SPH(1,:)/deg2rad; % theta in degrees
phi        = GCOORD_SPH(2,:)/deg2rad; % phi in degrees
r          = GCOORD_SPH(3,:);
points_ref = double(theta >= theta_ref_n  & theta <= theta_ref_s & ...
                    phi   >= phi_ref_w    & phi   <= phi_ref_e & ...
                    r     >= r_surf-d_ref)';
end % END OF SUBFUNCTION compute_els_inside_refined_region