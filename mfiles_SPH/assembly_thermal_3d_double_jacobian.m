function [CC,rhs] = assembly_thermal_3d_double_jacobian(MESH,EL2NOD,PhaseID,T,dt,VAR,num_scaling,theta_cone,PHYSICS,OPTS_T,fidl)
% Usage: [CC,rhs] = assembly_thermal_3d_double_jacobian(MESH,EL2NOD,PhaseID,T,dt,VAR,num_scaling,theta_cone,PHYSICS,OPTS_T,fidl)
% 
% Purpose: 
%   Calculate the element matrices and vectors; assemble global
%   counterparts 
%
% Input:
%   MESH        : [structure] : FE mesh parameters
%   EL2NOD      : [matrix]    : connectivity matrix
%   PhaseID     : [vector]    : phase indices
%   T           : [colvector] : variable field to be diffused
%   dt          : [scalar]    : time over which is diffused
%   VAR         : [structure] : major variable fields
%   num_scaling : [scalar]    : scaling factor to macth the non-SI units
%   theta_cone  : [scalar]    : half aperture of the cone
%   PHYSICS     : [structure] : physical properties
%   OPTS_D      : [structure] : options for assembly procedure
%   fidl        : [scalar]    : 
% 
% Output:
%   CC          : [matrix]    : global conductivity matrix
%   Rhs         : [vector]    : right hand side vector
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Aug 2013
% JMT Nov 2016 : double jacobian version (only 'std' method)

switch OPTS_T.method
    case 'std' % standard element assembly: good for understanding FE
        [CC,rhs] = assembly_std...
            (MESH,EL2NOD,PhaseID,T,dt,VAR,num_scaling,theta_cone,PHYSICS,OPTS_T,fidl); % *SUBFUNCTION*
    case 'opt' % a faster block-wise assembly (MILAMIN style)
        [CC,rhs] = assembly_opt...
            (MESH,EL2NOD,PhaseID,T,dt,VAR,num_scaling,theta_cone,PHYSICS,OPTS_T,fidl); % *SUBFUNCTION*
% %         [CC2,rhs2] = assembly_std...
% %             (MESH,EL2NOD,PhaseID,T,dt,VAR,num_scaling,theta_cone,PHYSICS,OPTS_T,fidl); % *SUBFUNCTION*
% %         disp(full(max(max(abs(CC-CC2))) / max(max((abs(CC))))))
% %         disp(max(abs(rhs-rhs2)) / max(abs(rhs)))
end

end % END OF FUNCTION assembly_thermal_3d_double_jacobian

% #########################################################################
%                               SUBFUNCTIONS
% #########################################################################

function [CC,rhs] = assembly_std(MESH,EL2NOD,PhaseID,T,dt,VAR,num_scaling,theta_cone,PHYSICS,OPTS_T,fidl)

t      = tic;

%==========================================================================
% MODEL INFO
%==========================================================================
nel    = size(EL2NOD,2);
nnodel = size(EL2NOD,1);
nnod   = length(T);

%==========================================================================
% COMPUTE ELEMENTS IN RELATION WITH THE CONE AND CROSSING PHI = 2PI
%==========================================================================
if strcmp(OPTS_T.eltype,'linear')
        % Compute them again since EL2NOD is now for 4-node linear elements
        % (after spliting each quadratic 10-node element into 8 linear
        % 4-node elements). The data stored in MESH is regarding a 10-node
        % element mesh
        [els_in_cone,~,~,els_out_cone,els_in_cone_iso] = ...
            check_els_in_cone(MESH.GCOORD_SPH,EL2NOD,theta_cone);
        [els_cross_2pi,~]         = check_phi(MESH.GCOORD_SPH,EL2NOD(:,els_out_cone));
        els_out_cone_cross_2pi    = els_out_cone(els_cross_2pi);
        els_out_cone_no_cross_2pi = els_out_cone;
        els_out_cone_no_cross_2pi(ismember(els_out_cone_no_cross_2pi,els_out_cone_cross_2pi)) = [];
        els_in_cone_no_iso        = els_in_cone;
        els_in_cone_no_iso(ismember(els_in_cone_no_iso,els_in_cone_iso)) = [];
        
        MESH.els_out_cone_no_cross_2pi = els_out_cone_no_cross_2pi;
        MESH.els_out_cone_cross_2pi    = els_out_cone_cross_2pi;
        MESH.els_in_cone_no_iso        = els_in_cone_no_iso;
        MESH.els_in_cone_iso           = els_in_cone_iso;
end

%==========================================================================
% PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
%==========================================================================
[x_ip,SF.w_ip] = ip_tetrahedron(OPTS_T.nip);
[SF.N,SF.dNds] = sf_dsf_tet(x_ip,nnodel,'cell');
[SF.N4,dN4ds]  = sf_dsf_tet(x_ip,4,'cell'); % linear s.f. and their derivatives
SF.dN4ds       = dN4ds{1}';                 % used to calculate each element's Jacobian

%==========================================================================
% STORAGE FOR GLOBAL MATRICES AND VECTORS
%==========================================================================
ASSEMBLY.CC_all = zeros(nnodel*nnodel,nel); % storage for conductivity matrix coefficients
ASSEMBLY.i_all  = zeros(nnodel*nnodel,nel); % storage for the matrix i indixes
ASSEMBLY.j_all  = zeros(nnodel*nnodel,nel); % storage for the matrix j indixes
ASSEMBLY.rhs    = zeros(nnod,1);            % storage for rhs force vector coefficients

%==========================================================================
% ASSEMBLY MATRICES FOR EACH SET OF ELEMENTS
%==========================================================================
ASSEMBLY = assembly_thermal_els_out_cone_no_cross_2pi_std...
    (MESH,EL2NOD,PhaseID,T,dt,VAR,num_scaling,ASSEMBLY,SF,PHYSICS,OPTS_T);
ASSEMBLY = assembly_thermal_els_out_cone_cross_2pi_std...
    (MESH,EL2NOD,PhaseID,T,dt,VAR,num_scaling,ASSEMBLY,SF,PHYSICS,OPTS_T);
ASSEMBLY = assembly_thermal_els_in_cone_no_iso_std...
    (MESH,EL2NOD,PhaseID,T,dt,VAR,num_scaling,ASSEMBLY,SF,PHYSICS,OPTS_T);
ASSEMBLY = assembly_thermal_els_in_cone_iso_std...
    (MESH,EL2NOD,PhaseID,T,dt,VAR,num_scaling,ASSEMBLY,SF,PHYSICS,OPTS_T);
t_elloop = toc(t);t=tic;

%==========================================================================
% MERGE ELEMENT MATRICES INTO GLOBAL SPARSE MATRIX
%==========================================================================
CC  = sparse3(ASSEMBLY.i_all(:),ASSEMBLY.j_all(:),ASSEMBLY.CC_all(:)); % global conductivity matrix
CC  = tril(CC); % this is to have the same CC than in 'opt' method
rhs = ASSEMBLY.rhs;
t_sparse = toc(t);

if fidl
    fprintf(fidl,'   element loop                          : %7.2f sec \n',t_elloop);
    fprintf(fidl,'   making sparse matrix                  : %7.2f sec\n',t_sparse);
end

end % END OF SUBFUNCTION assembly_std

% #########################################################################

function [CC,rhs] = assembly_opt(MESH,EL2NOD,PhaseID,T,dt,VAR,num_scaling,theta_cone,PHYSICS,OPTS_T,fidl)

% MILAMIN method (blocks of elements handled at once)
% see Dabrowski et al. 2008 (doi:10.1029/2007GC001719)

t      = tic;

%==========================================================================
% MODEL INFO
%==========================================================================
nel     = size(EL2NOD,2);
nnodel  = size(EL2NOD,1);

%==========================================================================
% COMPUTE ELEMENTS IN RELATION WITH THE CONE AND CROSSING PHI = 2PI
%==========================================================================
if strcmp(OPTS_T.eltype,'linear')
        % Compute them again since EL2NOD is now for 4-node linear elements
        % (after spliting each quadratic 10-node element into 8 linear
        % 4-node elements). The data stored in MESH is regarding a 10-node
        % element mesh
        [els_in_cone,~,~,els_out_cone,els_in_cone_iso] = ...
            check_els_in_cone(MESH.GCOORD_SPH,EL2NOD,theta_cone);
        [els_cross_2pi,~]         = check_phi(MESH.GCOORD_SPH,EL2NOD(:,els_out_cone));
        els_out_cone_cross_2pi    = els_out_cone(els_cross_2pi);
        els_out_cone_no_cross_2pi = els_out_cone;
        els_out_cone_no_cross_2pi(ismember(els_out_cone_no_cross_2pi,els_out_cone_cross_2pi)) = [];
        els_in_cone_no_iso        = els_in_cone;
        els_in_cone_no_iso(ismember(els_in_cone_no_iso,els_in_cone_iso)) = [];
        
        MESH.els_out_cone_no_cross_2pi = els_out_cone_no_cross_2pi;
        MESH.els_out_cone_cross_2pi    = els_out_cone_cross_2pi;
        MESH.els_in_cone_no_iso        = els_in_cone_no_iso;
        MESH.els_in_cone_iso           = els_in_cone_iso;
end

%==========================================================================
% PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
%==========================================================================
[x_ip,SF.w_ip] = ip_tetrahedron(OPTS_T.nip);
[SF.N,SF.dNds] = sf_dsf_tet(x_ip,nnodel,'cell');
[SF.N4,dN4ds]  = sf_dsf_tet(x_ip,4,'cell'); % linear s.f. and their derivatives
SF.dN4ds       = dN4ds{1};                  % used to calculate each element's Jacobian

%==========================================================================
% DECLARE VARIABLES (ALLOCATE MEMORY)
%==========================================================================
ASSEMBLY.CC_v  = zeros(nnodel*(nnodel+1)/2,nel);
ASSEMBLY.rhs_v = zeros(nnodel,nel);
ASSEMBLY.nblk  = 0;

%==========================================================================
% ASSEMBLY MATRICES FOR EACH SET OF ELEMENTS
%==========================================================================
ASSEMBLY = assembly_thermal_els_out_cone_no_cross_2pi_opt...
    (MESH,EL2NOD,PhaseID,T,dt,VAR,num_scaling,ASSEMBLY,SF,PHYSICS,OPTS_T);
ASSEMBLY = assembly_thermal_els_out_cone_cross_2pi_opt...
    (MESH,EL2NOD,PhaseID,T,dt,VAR,num_scaling,ASSEMBLY,SF,PHYSICS,OPTS_T);
ASSEMBLY = assembly_thermal_els_in_cone_no_iso_opt...
    (MESH,EL2NOD,PhaseID,T,dt,VAR,num_scaling,ASSEMBLY,SF,PHYSICS,OPTS_T);
ASSEMBLY = assembly_thermal_els_in_cone_iso_opt...
    (MESH,EL2NOD,PhaseID,T,dt,VAR,num_scaling,ASSEMBLY,SF,PHYSICS,OPTS_T);
t_elloop = toc(t);t=tic;

%==========================================================================
% ASSEMBLE R.H.S. VECTOR
%==========================================================================
rhs = accumarray(EL2NOD(:), ASSEMBLY.rhs_v(:));

%==========================================================================
% MERGE ELEMENT MATRICES INTO GLOBAL SPARSE MATRIX
%==========================================================================
opts_mutils.symmetric  = 1;
opts_mutils.n_node_dof = 1;
opts_mutils.nthreads   = 1; % >1 is very slow (????)
CC   = sparse_create(EL2NOD,ASSEMBLY.CC_v,opts_mutils);
nblk = ASSEMBLY.nblk;
t_sparse = toc(t);

if fidl
    fprintf(fidl,'   element loop (nelblk=%4i,nblk=%6i): %7.2f sec \n',...
        OPTS_T.nelblk,nblk,t_elloop);
    fprintf(fidl,'   making sparse matrix                  : %7.2f sec\n',t_sparse);
end

end % END OF SUBFUNCTION assembly_opt