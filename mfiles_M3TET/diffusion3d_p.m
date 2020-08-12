function [T,CC,MM] = diffusion3d_p(T,kappa,MESH,COMM,SETTINGS,TBC,dt,CC,MM)
% Usage: [T,CC,MM] = diffusion3d_p(T,kappa,MESH,COMM,SETTINGS,TBC,dt,CC,MM)
% 
% Purpose: Parallel finite element diffusion solver
%
% Input:
%   T        : [colvector] : variable field to be diffused
%   kappa    : [scalar]    : diffusivity (must be in correct units !!!)
%   MESH     : [structure] : FE mesh parameters
%   COMM     : [structure] : inter-subdomain communication data
%   SETTINGS : [structure] : model parameters
%   TBC      : [structure] : boundary conditions
%   dt       : [scalar]    : time step (must be in correct units !!!)
%  optional (if CC and MM are provided they won't be assembled)
%   CC       : [matrix]    : global conductivity matrix
%   MM       : [matrix]    : global capacity matrix
% 
% Output:
%   T        : [colvector] : diffused variable field
%   CC       : [matrix]    : global conductivity matrix
%   MM       : [matrix]    : global capacity matrix
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Apr 2011
% JH Aug 2013
% JH Nov 2016 : cleaned up, restructured
%

fidl   = SETTINGS.fid_log;
OPTS_D = SETTINGS.OPTS_D;

if SETTINGS.disp_profiling
    fprintf(fidl,' DIFFUSION \n');
else
    fprintf(fidl,' DIFFUSION '); fidl = 0;
end

if length(unique(kappa))>1
    error('Diffusivity kappa must be a scalar.');
end
if kappa==0 % No need to solve for diffusion if kappa==0
    if nargin<9
        CC = [];
        MM = [];
    end
    return
end


%==========================================================================
% MODEL INFO
%==========================================================================
if iscell(MESH.EL2NOD)
    EL2NOD = MESH.EL2NOD{1};
else
    EL2NOD = MESH.EL2NOD;
end
nnodel     = size(EL2NOD,1);
check_options(OPTS_D,nnodel); % *SUBFUNCTION*
if nnodel==10 && strcmp(OPTS_D.eltype,'linear')
    EL2NOD = tetmesh_p2_to_p1(MESH.GCOORD,EL2NOD);
end


%==========================================================================
% ASSEMBLY OF GLOBAL MATRICES (CC and MM contain BOTH SYMMETRIC HALVES)
%==========================================================================
switch SETTINGS.jacobian
    case 'standard'
        [CC,MM] = assembly_diffusion_3d(MESH.GCOORD,EL2NOD,kappa,OPTS_D,fidl);
    case 'double'
        [CC,MM] = assembly_diffusion_3d_double_jacobian(MESH,EL2NOD,kappa,SETTINGS.theta_cone,OPTS_D,fidl);
end

if strcmp(OPTS_D.lumping,'yes')
    % Lump mass matrix to improve stability
    MM = diag(sum(MM,2));
end


%==========================================================================
% Time approximation scheme (theta-time stepping rule)
%==========================================================================
% theta == 0   :: forward difference (EXPLICIT)
%       == 1   :: backward difference (FULL IMPLICIT)
%       == 1/2 :: midpoint rule; Crank-Nicholson (IMPLICIT)
%       == 2/3 :: Galerkin (IMPLICIT)
theta = OPTS_D.theta;
rhs   = MM*T - dt*(1-theta)* CC * T;
CC    = dt*theta*CC + MM;


%==========================================================================
% BOUNDARY CONDITIONS
%==========================================================================
nnod  = length(T);
ifree = setdiff(1:nnod,TBC.ifix); % replaces Free = 1:nnod; Free(TBC..nod)= [];
if ~isempty(TBC.ifix)
    rhs = rhs - CC(:,TBC.ifix)*TBC.vfix';
end
if isfield(TBC,'NN_INTEG')
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
T(ifree) = solve_matrix_eq(CC(ifree,ifree),rhs(ifree),T(ifree),COMM,OPTS_D,fidl);

end % END OF FUNCTION diffusion3d_p

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function check_options(OPTS_D,nnodel)

required_fields = {'solver' 'theta' 'eltype' 'lumping'};
for i=1:length(required_fields)
    if ~isfield(OPTS_D,required_fields{i})
        error('OPTS_D.%s is required in diffusion3d_p',required_fields{i});
    end
end

if strcmp(OPTS_D.eltype,'quadratic') 
    if nnodel==4
        error('Need 10-node tetrahedra mesh when using quadratic elements.');
    end
    if strcmp(OPTS_D.lumping,'yes')
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
    
end % END OF SUBFUNCTION include_BCs