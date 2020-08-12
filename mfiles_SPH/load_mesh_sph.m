function [GCOORD,EL2NOD,PointID,PhaseID,r_cmb,r_surf,EMBEDDED_MESH,el2sd] = ...
    load_mesh_sph(filename,nnodel,fidl,DB_indices)

tmp      = load(filename);
MESH     = tmp.MESH; clear tmp
GCOORD   = MESH.GCOORD;
EL2NOD   = uint32(MESH.EL2NOD);

if ~isfield(MESH,'r_cmb')
    MESH.r_cmb = 3471;
end
if ~isfield(MESH,'r_surf')
    MESH.r_surf = 6371;
end

% % % Check if we have nodes that are very close to each other (double nodes)
% % % (very slow for larger meshes, only use if necessary)
% % [D,I]     = pdist2(GCOORD(:,:)',GCOORD(:,:)','euclidean','SMALLEST',2);
% % el_length = calc_tetra_length(GCOORD,EL2NOD);
% % [mindist_nods,inod] = min(D(2,:));
% % if mindist_nods<0.1*min(el_length)
% %     error('Distance between nodes %1i and %1i: %.4e km',inod,I(2,inod),mindist_nods);
% % end
% % clear D I el_length

% Check if node connectivity is ok; If a negative volume is calculated
% (when assuming Hughes' notation), convert the numbering to Hughes' notation 
el_vol        = calc_tetra_volume(GCOORD,EL2NOD);
iel           = el_vol<0;
if any(iel)
    EL2NOD(1:4,iel) = MESH.EL2NOD([1 2 4 3],iel); % change connectivity
end

% Load PointID
PointID = int32(MESH.PointID);

% Load PhaseID (or set to uniform "1" everywhere)
if isfield(MESH,'PhaseID')
    PhaseID = int32(MESH.PhaseID);
else
    PhaseID = ones(1,size(EL2NOD,2),'int32');
end

% Load radius for CMB and surface
r_cmb  = MESH.r_cmb;
r_surf = MESH.r_surf;

try MESH.l0_ref; % the loaded mesh has an embedded high resolution subregion
    EMBEDDED_MESH.l0_ref = MESH.l0_ref; % element length inside refined region
    EMBEDDED_MESH.theta0 = MESH.theta0; % colatitude (degrees) of the point around which the refined zone is defined
    EMBEDDED_MESH.phi0   = MESH.phi0;   % longitude (degrees) of the point around which the refined zone is defined
    EMBEDDED_MESH.d_ref  = MESH.d_ref;  % depth
    EMBEDDED_MESH.w_ref  = MESH.w_ref;  % width (North-South)
    EMBEDDED_MESH.l_ref  = MESH.l_ref;  % length (East-West)
    EMBEDDED_MESH.d_tran = MESH.d_tran; % depth
    EMBEDDED_MESH.w_tran = MESH.w_tran; % width (North-South)
    EMBEDDED_MESH.l_tran = MESH.l_tran; % length (East-West)
    % Finite rotation matrices to change the reference frame between GPlates and the spherical mesh
    EMBEDDED_MESH.RR1    = MESH.RR1_center;
    EMBEDDED_MESH.RR2    = MESH.RR2_center;
catch
    EMBEDDED_MESH.l0_ref = [];
    EMBEDDED_MESH.theta0 = [];
    EMBEDDED_MESH.phi0   = [];
    EMBEDDED_MESH.d_ref  = [];
    EMBEDDED_MESH.w_ref  = [];
    EMBEDDED_MESH.l_ref  = [];
    EMBEDDED_MESH.d_tran = [];
    EMBEDDED_MESH.w_tran = [];
    EMBEDDED_MESH.l_tran = [];
    EMBEDDED_MESH.RR1    = [];
    EMBEDDED_MESH.RR2    = [];
end

% Check if a suggested subdomain split is provided ("el2sd" is a pointer element -->subdomain)
if numlabs>1 && isfield(MESH,['el2sd_' num2str_d(numlabs,3)])
    el2sd = MESH.(['el2sd_' num2str_d(numlabs,3)]);
else
    el2sd = [];
end

% Add nodes or re-connect nodes if necessary
switch nnodel
    case 4
        if size(EL2NOD,1)==10
            % Create 4-node connectivity matrix by reconnecting the nodes
            fprintf(fidl,'\n Quadratic (10-node) connectivity was loaded from file. \n');
            [EL2NOD,PhaseID] = tetmesh_p2_to_p1(GCOORD,EL2NOD,PhaseID);
            fprintf(fidl,' Linear (4-node) connectivity has been generated. \n');
        end
        
    case 10
        if size(EL2NOD,1)==4
            % Create edge nodes (10-node elements) if a linear (4-node) mesh was loaded
            fprintf(fidl,'\n Linear (4-node) connectivity was loaded from file. \n');
            [GCOORD,EL2NOD,PointID] = tetmesh_p1_to_p2(GCOORD,EL2NOD,PointID,DB_indices);
            fprintf(fidl,' Quadratic (10-node) connectivity has been generated. \n');
        end
        
    otherwise
        error(' nnodel must be either 4 or 10.');
end

end % END OF FUNCTION load_mesh_sph