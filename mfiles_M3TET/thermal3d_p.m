function T = thermal3d_p(T,VAR,MESH,COMM,PHYSICS,SETTINGS,NUMSCALE,TBC,dt)
% Usage: T = thermal3d_p(T,VAR,MESH,COMM,PHYSICS,SETTINGS,NUMSCALE,TBC,dt)
% 
% Purpose: Parallel fully implicit finite element thermal diffusion solver
%
% Input:
%   T        : [colvector] : variable field to be diffused
%   VAR      : [structure] : major variable fields
%   MESH     : [structure] : FE mesh parameters
%   COMM     : [structure] : inter-subdomain communication data
%   PHYSICS  : [structure] : physical properties
%   SETTINGS : [structure] : model parameters
%   NUMSCALE : [structure] : numerical scaling parameters
%   TBC      : [structure] : boundary conditions
%   dt       : [scalar]    : time over which is diffused
% 
% Output:
%   T        : [colvector] : diffused variable field
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Aug 2013
% JH Nov 2016 : cleaned up, restructured
%

fidl   = SETTINGS.fid_log;
OPTS_T = SETTINGS.OPTS_T;

if SETTINGS.disp_profiling
    fprintf(fidl,' TEMPERATURE \n');
else
    fprintf(fidl,' TEMPERATURE '); fidl = 0;
end


%==========================================================================
% MODEL INFO
%==========================================================================
if iscell(MESH.EL2NOD)
    EL2NOD = MESH.EL2NOD{1};
else
    EL2NOD = MESH.EL2NOD;
end
PhaseID    = MESH.PhaseID;
nnodel     = size(EL2NOD,1);
check_options(OPTS_T,nnodel); % *SUBFUNCTION*
if nnodel==10 && strcmp(OPTS_T.eltype,'linear')
    [EL2NOD,PhaseID] = tetmesh_p2_to_p1(EL2NOD,PhaseID);
end


%==========================================================================
% USE SCALING FACTOR TO MATCH THE NON-SI UNITS
%==========================================================================
num_scaling = (NUMSCALE.t0/NUMSCALE.L0^2); % scaling factor
% if num_scaling~=1
%     % CONVERT NON-SI-UNIT VARIABLES TO SI-UNITS
%     GCOORD = GCOORD * NUMSCALE.L0; % km --> m
%     dt     = dt * NUMSCALE.t0;     % Myr --> s
%     num_scaling = 1;
% end


%==========================================================================
% ASSEMBLY OF GLOBAL MATRICES
%==========================================================================
switch SETTINGS.jacobian
    case 'standard'
        [CC,rhs] = assembly_thermal_3d(MESH.GCOORD,EL2NOD,PhaseID,T,dt,VAR,num_scaling,PHYSICS,OPTS_T,fidl);
    case 'double'
        [CC,rhs] = assembly_thermal_3d_double_jacobian(MESH,EL2NOD,PhaseID,T,dt,VAR,num_scaling,SETTINGS.theta_cone,PHYSICS,OPTS_T,fidl);
end


%==========================================================================
% BOUNDARY CONDITIONS
%==========================================================================
nnod  = length(T);
ifree = setdiff(1:nnod,TBC.ifix); % replaces Free = 1:nnod; Free(TBC..nod)= [];
if ~isempty(TBC.ifix)
    TMP = CC(:,TBC.ifix) + cs_transpose(CC(TBC.ifix,:)); % diagonal values twice ?????
    rhs = rhs - TMP*TBC.vfix'; clear TMP
end
if isfield(TBC,'NN_INTEG')
    if num_scaling~=1
        error('Heat flux boundary conditions must be verified if num_scaling~=1 !!');
    end
    valHF           = zeros(nnod,1);
    valHF(TBC.iflx) = TBC.vflx;
    HF              = num_scaling*dt*TBC.NN_INTEG*valHF; % FIXME now scaling with num_scaling to make similar to thermal2d
    rhs             = rhs + HF; % add heat flux to rhs
end
T(TBC.ifix) = TBC.vfix;


% =========================================================================
% Update communication arrays to account for the new size of the system
% of equations to be solved (ifree instead of nnod).
% =========================================================================
if COMM.nsd>1
    COMM = include_BCs(COMM,TBC.ifix,ifree);
end


%==========================================================================
% Solve equation
%==========================================================================
if SETTINGS.disp_profiling
    fprintf(fidl,'   solving eq   ');
end
T(ifree) = solve_matrix_eq(CC(ifree,ifree),rhs(ifree),T(ifree),COMM,OPTS_T,fidl);

end % END OF FUNCTION thermal3d_p

% #########################################################################
%                               SUBFUNCTIONS
% #########################################################################

function check_options(OPTS_T,nnodel)

required_fields = {'method_eval_dens' 'method_eval_cond' ...
                   'method_eval_Cp' 'method_eval_dQdt' ...
                   'solver' 'eltype' 'lumping'};
for i=1:length(required_fields)
    if ~isfield(OPTS_T,required_fields{i})
        error('OPTS_T.%s is required in thermal3d_p',required_fields{i});
    end
end

if strcmp(OPTS_T.eltype,'quadratic') 
    if nnodel==4
        error('Need 10-node tetrahedra mesh when using quadratic elements.');
    end
    if strcmp(OPTS_T.lumping,'yes')
        error('Lumping only works with 4-node tetrahedra.');
    end
end

end % END OF SUBFUNCTION check_options

% #########################################################################

function COMM = include_BCs(COMM,ifix,ifree)

% (1) construct pointer from free dofs to all (free+fixed) dofs
nfree                = length(ifree);
free2all_dofs(ifree) = 1:nfree;

% (2) Use the pointer to extract and re-number the free dofs
nNB              = length(COMM.NB);
mynodSDB         = COMM.mynod_SDB{1};
COMM.mynod_SDB_f = cell(1,nNB);
for iNB=1:nNB
    if COMM.NB(iNB)
        % Extract the SD boundary dofs shared with SD "NB(iNB)" but without
        % changing their order!
        mynodSDB_free = mynodSDB{iNB}( ~ismember(mynodSDB{iNB},ifix) );
        % Use the above pointer to re-number the SD boundary dofs
        % The dofs now correspond to the reduced system of equations:
        % KK(free,free) * u(free) = rhs(free)
        COMM.mynod_SDB_f{iNB} = free2all_dofs( mynodSDB_free );
    end
end

% (3) Same re-numbering for the list of uniq dofs on finest MG level
unique_nodes_f = COMM.unique_nodes{1}( ~ismember(COMM.unique_nodes{1},ifix) );
COMM.unique_nodes_f = free2all_dofs( unique_nodes_f );

% (4) Calculate the number of unknown temperatures in whole domain
%     (needed for norm calculation)
COMM.nD = COMM.sum_all( length(unique_nodes_f) );
    
end