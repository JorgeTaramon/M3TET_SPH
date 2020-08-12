function [CC,MM] = assembly_diffusion_3d_double_jacobian(MESH,EL2NOD,kappa,theta_cone,OPTS_D,fidl)
% Usage: [CC,MM] = assembly_diffusion_3d_double_jacobian(MESH,EL2NOD,kappa,theta_cone,OPTS_D,fidl)
%
% Purpose: 
%   Calculate the element matrices and vectors; assemble global
%   counterparts 
%
% Input:
%   MESH       : [structure] : FE mesh parameters
%   EL2NOD     : [matrix]    : connectivity matrix
%   kappa      : [scalar]    : diffusivity (must be in correct units !!!)
%   theta_cone : [scalar]    : half aperture of the cone
%   OPTS_D     : [structure] : options for assembly procedure
%   fidl       : [scalar]    : 
%
% Output:
%   CC         : [sparsemat] : global conductivity (diffusivity) matrix
%   MM         : [sparsemat] : global mass (capacity) matrix
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)

% JH March 2011
% JH March 2013
% JMT Nov 2016 : double jacobian version (only 'std' method)

switch OPTS_D.method
    case 'std' % standard element assembly: good for understanding FE
        [CC,MM] = assembly_std...
            (MESH,EL2NOD,kappa,theta_cone,OPTS_D,fidl); % *SUBFUNCTION*
    case 'opt' % a faster block-wise assembly (MILAMIN style)
        [CC,MM] = assembly_opt...
            (MESH,EL2NOD,kappa,theta_cone,OPTS_D,fidl); % *SUBFUNCTION*
% %         [CC2,MM2] = assembly_std...
% %             (MESH,EL2NOD,kappa,theta_cone,OPTS_D,fidl); % *SUBFUNCTION*
% %         disp(full(max(max(abs(CC-CC2))) / max(max((abs(CC))))))
% %         disp(full(max(max(abs(MM-MM2))) / max(max((abs(MM))))))
end

end % END OF FUNCTION assembly_diffusion_3d_double_jacobian

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function [CC,MM] = assembly_std(MESH,EL2NOD,kappa,theta_cone,OPTS_D,fidl)

t      = tic;

%==========================================================================
% MODEL INFO
%==========================================================================
nel        = size(EL2NOD,2);
nnodel     = size(EL2NOD,1);

%==========================================================================
% COMPUTE ELEMENTS IN RELATION WITH THE CONE AND CROSSING PHI = 2PI
%==========================================================================
if strcmp(OPTS_D.eltype,'linear')
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
[x_ip,SF.w_ip] = ip_tetrahedron(OPTS_D.nip);
[SF.N,SF.dNds] = sf_dsf_tet(x_ip,nnodel,'cell');
[SF.N4,dN4ds]  = sf_dsf_tet(x_ip,4,'cell'); % linear s.f. and their derivatives
SF.dN4ds       = dN4ds{1}';                    % used to calculate each element's Jacobian

%==========================================================================
% STORAGE FOR GLOBAL MATRICES AND VECTORS
%==========================================================================
ASSEMBLY.CC_all = zeros(nnodel*nnodel,nel); % storage for conductivity matrix coefficients
ASSEMBLY.MM_all = zeros(nnodel*nnodel,nel); % storage for capacity matrix coefficients
ASSEMBLY.i_all  = zeros(nnodel*nnodel,nel); % storage for the matrix i indixes
ASSEMBLY.j_all  = zeros(nnodel*nnodel,nel); % storage for the matrix j indixes

%==========================================================================
% ASSEMBLY MATRICES FOR EACH SET OF ELEMENTS
%==========================================================================
ASSEMBLY = assembly_diffusion_els_out_cone_no_cross_2pi_std...
    (MESH,EL2NOD,ASSEMBLY,kappa,SF,OPTS_D);
ASSEMBLY = assembly_diffusion_els_out_cone_cross_2pi_std...
    (MESH,EL2NOD,ASSEMBLY,kappa,SF,OPTS_D);
ASSEMBLY = assembly_diffusion_els_in_cone_no_iso_std...
    (MESH,EL2NOD,ASSEMBLY,kappa,SF,OPTS_D);
ASSEMBLY = assembly_diffusion_els_in_cone_iso_std...
    (MESH,EL2NOD,ASSEMBLY,kappa,SF,OPTS_D);
t_elloop = toc(t);t=tic;

%==========================================================================
% MERGE ELEMENT MATRICES INTO GLOBAL SPARSE MATRIX
%==========================================================================
CC  = sparse3(ASSEMBLY.i_all(:),ASSEMBLY.j_all(:),ASSEMBLY.CC_all(:)); % global conductivity matrix
MM  = sparse3(ASSEMBLY.i_all(:),ASSEMBLY.j_all(:),ASSEMBLY.MM_all(:)); % global conductivity matrix
t_sparse = toc(t);

if fidl
    fprintf(fidl,'   element loop                          : %7.2f sec \n',t_elloop);
    fprintf(fidl,'   making sparse matrix                  : %7.2f sec\n',t_sparse);
end

end % END OF SUBFUNCTION assembly_std

% #########################################################################

function [CC,MM] = assembly_opt(MESH,EL2NOD,kappa,theta_cone,OPTS_D,fidl)

% MILAMIN method (blocks of elements handled at once)
% see Dabrowski et al. 2008 (doi:10.1029/2007GC001719)

t      = tic;

%==========================================================================
% MODEL INFO
%==========================================================================
nel     = size(EL2NOD,2); % number of elements
nnodel  = size(EL2NOD,1); % number of nodes in each element

%==========================================================================
% COMPUTE ELEMENTS IN RELATION WITH THE CONE AND CROSSING PHI = 2PI
%==========================================================================
if strcmp(OPTS_D.eltype,'linear')
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
[x_ip,SF.w_ip] = ip_tetrahedron(OPTS_D.nip);
[SF.N,SF.dNds] = sf_dsf_tet(x_ip,nnodel,'cell');
[SF.N4,dN4ds]  = sf_dsf_tet(x_ip,4,'cell'); % linear s.f. and their derivatives
SF.dN4ds       = dN4ds{1};                     % used to calculate each element's Jacobian

%==========================================================================
% DECLARE VARIABLES (ALLOCATE MEMORY)
%==========================================================================
ASSEMBLY.CC_v  = zeros(nnodel*(nnodel+1)/2,nel); % storage for conductivity matrix coefficients
ASSEMBLY.MM_v  = zeros(nnodel*(nnodel+1)/2,nel); % storage for capacity matrix coefficients
ASSEMBLY.nblk  = 0;

%==========================================================================
% ASSEMBLY STIFFNESS MATRICES FOR EACH SET OF ELEMENTS
%==========================================================================
ASSEMBLY = assembly_diffusion_els_out_cone_no_cross_2pi_opt...
    (MESH,EL2NOD,ASSEMBLY,kappa,SF,OPTS_D);
ASSEMBLY = assembly_diffusion_els_out_cone_cross_2pi_opt...
    (MESH,EL2NOD,ASSEMBLY,kappa,SF,OPTS_D);
ASSEMBLY = assembly_diffusion_els_in_cone_no_iso_opt...
    (MESH,EL2NOD,ASSEMBLY,kappa,SF,OPTS_D);
ASSEMBLY = assembly_diffusion_els_in_cone_iso_opt...
    (MESH,EL2NOD,ASSEMBLY,kappa,SF,OPTS_D);
t_elloop = toc(t);t=tic;

%==========================================================================
% MERGE ELEMENT MATRICES INTO GLOBAL SPARSE MATRIX
%==========================================================================
opts_mutils.symmetric  = 1;
opts_mutils.n_node_dof = 1;
opts_mutils.nthreads   = 1; % >1 is very slow (????)
CC   = sparse_create(EL2NOD,ASSEMBLY.CC_v,opts_mutils);
MM   = sparse_create(EL2NOD,ASSEMBLY.MM_v,opts_mutils);
nblk = ASSEMBLY.nblk;
clear ASSEMBLY
t_sparse = toc(t);

% Create second half symmetric matrices
CC = CC + tril(CC,-1)';
MM = MM + tril(MM,-1)';

if fidl
    fprintf(fidl,'   element loop (nelblk=%4i,nblk=%6i): %7.2f sec \n',...
        OPTS_D.nelblk,nblk,t_elloop);
    fprintf(fidl,'   making sparse matrix                  : %7.2f sec\n',t_sparse);
end

end % END OF SUBFUNCTION assembly_opt