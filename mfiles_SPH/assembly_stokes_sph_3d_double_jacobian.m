function [KK,Fb,GG,MM,Fpdil,DensEl,ViscIP] = assembly_stokes_sph_3d_double_jacobian...
             (VAR,MESH,PHYSICS,NUMSCALE,OPTS)
% Usage: [KK,Fb,GG,MM,Fpdil,DensEl,ViscIP] = assembly_stokes_sph_3d_double_jacobian...
%            (VAR,MESH,PHYSICS,NUMSCALE,OPTS)
% 
% Purpose: Calculate element matrices and vectors for viscous flow problem 
%          and assemble global counterparts.
%
% Input:
%   VAR      : [structure] : major variable fields, each is a vector
%   MESH     : [structure] : FE mesh parameters
%   PHYSICS  : [structure] : physical properties
%   NUMSCALE : [structure] : numerical scaling parameters
%   OPTS     : [structure] : options for assembly procedure
%
% Output:
%   KK     : [sparsemat] : global stiffness matrix
%   Fb     : [colvector] : global (buoyancy) force vector
%   GG     : [sparsemat] : global gradient matrix
%   MM     : [sparsemat] : (1/viscosity)-scaled pressure mass matrix
%   Fpdil  : [colvector] : global dilation-vector (rhs of incompressibility
%                          equations)
%   DensEl : [matrix]    : density of all elements
%   ViscIP : [matrix]    : viscosity at integration points of all elements
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Jan 2011
% JPM 9Feb2011
% JH Jan 2012 : added switch for viscosity interp in elements
% JH Apr 2014 : cleaned up and made similar to visco-elastic assembly
%               function 
% JMT/JPM Mar 2016: matrices are computed for curved edges (in Cartesian
%                   coord) elements using a double mapping from local to
%                   spherical coordinates and then from spherical to
%                   Cartesian coordinates (only 'std' method). The logic
%                   takes into account that midside nodes of some elements
%                   (elements inside the cone) were computed in a spherical
%                   rotated frame (meaning straight-edge elements in that
%                   rotated frame).
% JH Nov 2016 : 'opt' version

% =========================================================================
% CHECK OPTIONS AND USE DEFAULT VALUES IF NECESSARY
% =========================================================================
if nargin<5
    error(' No input "OPTS" with assembly options defined.');
end
OPTS = check_options(OPTS); % *SUBFUNCTION*

% =========================================================================
% MATRIX ASSEMBLY
% =========================================================================
switch OPTS.method
    case 'std' % standard element assembly: good for understanding FE
        [KK,Fb,GG,MM,Fpdil,DensEl,ViscIP] = element_assembly_std...
            (VAR,MESH,PHYSICS,NUMSCALE,OPTS); % *SUBFUNCTION*
    case 'opt' % a faster block-wise assembly (MILAMIN style)
        [KK,Fb,GG,MM,Fpdil,DensEl,ViscIP] = element_assembly_opt...
            (VAR,MESH,PHYSICS,NUMSCALE,OPTS); % *SUBFUNCTION*
% %         [KK2,Fb2,GG2,MM2,Fpdil2,DensEl2,ViscIP2] = element_assembly_std...
% %             (VAR,MESH,PHYSICS,NUMSCALE,OPTS); % *SUBFUNCTION*
% %         max(max(abs(KK-KK2)))./max(max(abs(KK2)))
% %         max(max(abs(MM-MM2)))./max(max(abs(MM2)))
% %         max(max(abs(GG-GG2)))./max(max(abs(GG2)))
% %         max(abs(Fpdil-Fpdil2))./max(abs(Fpdil2))
% %         max(abs(Fb-Fb2))./max(abs(Fb2))
end

end % END OF FUNCTION assembly_stokes_sph_3d_double_jacobian

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function OPTS = check_options(OPTS)

% Method for assembly
if isfield(OPTS,'method')
    if ~strcmp(OPTS.method,'std') && ~strcmp(OPTS.method,'opt')
        error(' OPTS.method must be either "std" or "opt".');
    end
else
    OPTS.method = 'opt'; % DEFAULT (optimized assembly)
end

% Blocksize if optimized assembly
if strcmp(OPTS.method,'opt')
    if ~isfield(OPTS,'nelblk')
        OPTS.nelblk = 1000;% DEFAULT (number of elements per block)
    end
end

% Number of blocks that are assembled at once using sparse command
if ~isfield(OPTS,'nblk_asmbl') || isempty(OPTS.nblk_asmbl)
    % Check if SuiteSparse is installed so that the faster sparse2 can be
    % used instead of MATLAB's sparse-function
    try
        sparse2(1,1,1);
        OPTS.nblk_asmbl = ceil(40000/OPTS.nelblk);
    catch %#ok<CTCH>
        OPTS.nblk_asmbl = ceil(10000/OPTS.nelblk);
    end
end

% Threshold for values to be assembled
% (values smaller than rtol_assembly*max(abs(KKv)) won't be assembled)
if isfield(OPTS,'rtol_asmbly') && ~isempty(OPTS.rtol_asmbly)
    if OPTS.rtol_asmbly>1e-16
        error(' Threshold for values to be assembled seems to be too high');
    end
else
    OPTS.rtol_asmbly = 0;
end

% Use MUTILS' "sparse_create" instead of sparse2 or sparse
if ~isfield(OPTS,'use_mutils')
    OPTS.use_mutils = 0;
end

if OPTS.use_mutils
    if ~isfield(OPTS,'nthreads')
        OPTS.nthreads = 1;
    end
end

% Viscosity in each element
if ~isfield(OPTS,'method_eval_visc')
    error('OPTS.method_eval_visc must be defined.');
end

% Viscosity in each element
if ~isfield(OPTS,'method_eval_dens')
    error('OPTS.method_eval_dens must be defined.');
end

% D-matrix formulation
if ~isfield(OPTS,'type_Dmat')
    OPTS.type_Dmat = '221'; % '221' or '23rd'
end

end % END OF SUBFUNCTION check_options

% #########################################################################

function [KK,Fb,GG,MM,Fpdil,DensEl,ViscIP] = element_assembly_std...
             (VAR,MESH,PHYSICS,NUMSCALE,OPTS)

% Standard method: loop over integration points within loop over elements
% (becomes very slow for larger numbers of elements)

% =========================================================================
% MODEL PARAMETERS
% =========================================================================
EL2NOD      = double(MESH.EL2NOD{1});
EL2NODP     = double(MESH.EL2NOD{1}(1:4,:));
nnod        = MESH.nnod;
nel         = MESH.nel;
ndim        = 3;
nUdof       = nnod * ndim;
nPdof       = MESH.nVnod; % = max(max(EL2NODP))
nnodel      = size(EL2NOD,1);
nUdofel     = ndim*nnodel;
nPdofel     = size(EL2NODP,1);


%==========================================================================
% PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
%==========================================================================
[IP_X,SF.IP_w]       = ip_tetrahedron(OPTS.nip);
    % local coordinates and weights of points for integration of velocity/pressure matrices
[SF.NU,SF.dNUds]     = sf_dsf_tet(IP_X,nnodel,'cell');
    % velocity shape functions and their derivatives
[SF.NP,~]            = sf_dsf_tet(IP_X,nPdofel,'cell');
    % linear pressure shape functions and their derivatives
[~,dN4ds]            = sf_dsf_tet(IP_X,nPdofel,'cell');
SF.dN4ds             = dN4ds{1}';
    % derivatives of linear (4-node) shape functions (used to calculate each element's Jacobian 
    % when the elements have straight edges). These are the same for any global coordinate system
[SF.NU14,SF.dNUds14] = sf_dsf_tet(IP_X,size(MESH.EL2NOD14,1),'cell');
    % shape functions and their derivatives for 14-node elements. They are 
    % used to calculate the Jacobian for isoparametric elements
[SF.NU20,SF.dNUds20] = sf_dsf_tet(IP_X,size(MESH.EL2NOD_cubic,1),'cell');
    % shape functions and their derivatives for 20-node elements. They are 
    % used to calculate the Jacobian for isoparametric elements

%==========================================================================
% MATRIX FOR DEVIATORIC STRAIN RATE EXTRACTION
%==========================================================================
D = zeros(6);
switch OPTS.type_Dmat
    case '23rd'
        D(1,1) =  4/3; D(1,2) = -2/3; D(1,3) = -2/3; % See Zienkiewicz book, Vol 2
        D(2,1) = -2/3; D(2,2) =  4/3; D(2,3) = -2/3; % 4th edition, p. 519
        D(3,1) = -2/3; D(3,2) = -2/3; D(3,3) =  4/3; %
        D(4,4) =    1; D(5,5) =    1; D(6,6) =    1; %
    case '221'
        D(1,1) = 2; D(2,2) = 2; D(3,3) = 2; % See Zienkiewicz book, Vol 2
        D(4,4) = 1; D(5,5) = 1; D(6,6) = 1; % 4th edition, p. 335
end
B = zeros(6,nUdofel); % See Hughes' book p. 87

%==========================================================================
% STORAGE FOR GLOBAL MATRICES AND VECTORS
%==========================================================================
ASSEMBLY.Fb    = zeros(nUdof,1);             % storage for global force vector
ASSEMBLY.Fpdil = zeros(nPdof,1);             % storage for global dilation vector (pressure eq.)
ASSEMBLY.KKv   = zeros(nUdofel*nUdofel,nel); % storage for stiffness matrix data
ASSEMBLY.KKi   = zeros(nUdofel*nUdofel,nel); % storage for stiffness matrix i indices
ASSEMBLY.KKj   = zeros(nUdofel*nUdofel,nel); % storage for stiffness matrix j indices
ASSEMBLY.GGv   = zeros(nPdofel*nUdofel,nel); % storage for gradient matrix data
ASSEMBLY.GGi   = zeros(nPdofel*nUdofel,nel); % storage for gradient matrix i indices
ASSEMBLY.GGj   = zeros(nPdofel*nUdofel,nel); % storage for gradient matrix j indices

% Inverse-(viscosity*penalty)-scaled Mass matrix for pressure preconditioning
ASSEMBLY.MMv   = zeros(nPdofel*nPdofel,nel); % storage for (1/mu)-scaled mass matrix data
ASSEMBLY.MMi   = zeros(nPdofel*nPdofel,nel); % storage for mass matrix i indices
ASSEMBLY.MMj   = zeros(nPdofel*nPdofel,nel); % storage for mass matrix j indices

% Density and viscosity at integration points of all elements
if OPTS.return_dens_el
    ASSEMBLY.DensEl = zeros(nel,OPTS.nip);
else
    ASSEMBLY.DensEl = [];
end
if OPTS.return_visc_ip
    ASSEMBLY.ViscIP = zeros(nel,OPTS.nip);
else
    ASSEMBLY.ViscIP = [];
end

%==========================================================================
% ELEMENT TO GLOBAL VELOCITY D.O.F. CONNECTIVITY
%==========================================================================
EL2DOF                   = zeros(nUdofel, nel); % element -> global velocity dofs
EL2DOF(1:ndim:nUdofel,:) = ndim*EL2NOD-2;
EL2DOF(2:ndim:nUdofel,:) = ndim*EL2NOD-1;
EL2DOF(3:ndim:nUdofel,:) = ndim*EL2NOD;

%==========================================================================
% ASSEMBLY STIFFNESS MATRICES FOR EACH SET OF ELEMENTS
%==========================================================================
ASSEMBLY = assembly_stokes_els_out_cone_no_cross_2pi_std...
    (ASSEMBLY,SF,D,B,EL2DOF,VAR,MESH,PHYSICS,NUMSCALE,OPTS);
ASSEMBLY = assembly_stokes_els_out_cone_cross_2pi_std...
    (ASSEMBLY,SF,D,B,EL2DOF,VAR,MESH,PHYSICS,NUMSCALE,OPTS);
ASSEMBLY = assembly_stokes_els_in_cone_no_iso_std...
    (ASSEMBLY,SF,D,B,EL2DOF,VAR,MESH,PHYSICS,NUMSCALE,OPTS);
% ASSEMBLY = assembly_stokes_els_in_cone_iso_face_nodes_std...
%     (ASSEMBLY,SF,D,B,EL2DOF,VAR,MESH,PHYSICS,NUMSCALE,OPTS);
ASSEMBLY = assembly_stokes_els_in_cone_iso_20_nodes_std...
            (ASSEMBLY,SF,D,B,EL2DOF,VAR,MESH,PHYSICS,NUMSCALE,OPTS);

% %--------------------------------------------------------------------------
% % Uncomment following lines to use cubic isoparametric elements in the 
% % the rotated spherical frame (worse results!! -> IT NEEDS TO BE TESTED 
% % AS IT IS EXPECTED BETTER RESULTS) 
% % NOTE: Following functions need to be commented:
% %        - assembly_stokes_els_in_cone_iso_std
% %        - assembly_stokes_els_in_cone_iso_no_cross_2pi_std
% %        - assembly_stokes_els_in_cone_iso_cross_2pi_std
% %        - assembly_stokes_els_in_cone_iso_dealt_as_sub_std 
% % ASSEMBLY = assembly_stokes_cubic_els_in_cone_iso_std...
% %     (ASSEMBLY,D,B,EL2DOF,VAR,MESH,PHYSICS,NUMSCALE,OPTS);

% %--------------------------------------------------------------------------
% % Uncomment following lines to use isoparametric elements in the non
% % rotated spherical frame (litle worse results than using isoparametric 
% % elements in the rotated spherical frame)
% % NOTE: Following functions need to be commented:
% %        - assembly_stokes_els_in_cone_iso_std and
% %        - assembly_stokes_cubic_els_in_cone_iso_std
% %        - assembly_stokes_els_in_cone_iso_dealt_as_sub_std 
% [els_cross_2pi,~] = check_phi(MESH.GCOORD_SPH,MESH.EL2NOD{1}(1:4,MESH.els_in_cone_iso));
% MESH.els_in_cone_iso_cross_2pi    = MESH.els_in_cone_iso(els_cross_2pi);
% MESH.els_in_cone_iso_no_cross_2pi = MESH.els_in_cone_iso;
% MESH.els_in_cone_iso_no_cross_2pi(els_cross_2pi) = [];
% ASSEMBLY = assembly_stokes_els_in_cone_iso_no_cross_2pi_std...
%     (ASSEMBLY,SF,D,B,EL2DOF,VAR,MESH,PHYSICS,NUMSCALE,OPTS);
% ASSEMBLY = assembly_stokes_els_in_cone_iso_cross_2pi_std...
%     (ASSEMBLY,SF,D,B,EL2DOF,VAR,MESH,PHYSICS,NUMSCALE,OPTS);

% %--------------------------------------------------------------------------
% % Uncomment following lines to use isoparametric elements as if they were
% % subparametric in the rotated spherical frame (worse results!!)
% % NOTE: Following functions need to be commented:
% %        - assembly_stokes_els_in_cone_iso_std
% %        - assembly_stokes_els_in_cone_iso_no_cross_2pi_std
% %        - assembly_stokes_els_in_cone_iso_cross_2pi_std
% ASSEMBLY = assembly_stokes_els_in_cone_iso_dealt_as_sub_std...
%     (ASSEMBLY,SF,D,B,EL2DOF,VAR,MESH,PHYSICS,NUMSCALE,OPTS);

%==========================================================================
% ASSEMBLY OF GLOBAL SPARSE MATRICES AND RHS-VECTOR
%==========================================================================
MM = sparse3(ASSEMBLY.MMi(:),ASSEMBLY.MMj(:),ASSEMBLY.MMv(:)); % global (1/viscosity)-scaled mass matrix
GG = sparse3(ASSEMBLY.GGi(:),ASSEMBLY.GGj(:),ASSEMBLY.GGv(:)); % global gradient matrix
KK = sparse3(ASSEMBLY.KKi(:),ASSEMBLY.KKj(:),ASSEMBLY.KKv(:)); % global stiffness matrix
clear GGi GGj GGv MMi MMj MMv KKi KKj KKv

% only keep lower half for compatibility with optimized version
MM = tril(MM);
KK = tril(KK);

Fb     = ASSEMBLY.Fb;
Fpdil  = ASSEMBLY.Fpdil;
DensEl = ASSEMBLY.DensEl;
ViscIP = ASSEMBLY.ViscIP;

end % END OF FUNCTION element_assembly_std

% ###################################################################

function [KK,Fb,GG,MM,Fpdil,DensEl,ViscIP] = element_assembly_opt...
             (VAR,MESH,PHYSICS,NUMSCALE,OPTS)

% MILAMIN method (blocks of elements handled at once)
% see Dabrowski et al. 2008 (doi:10.1029/2007GC001719)

% =========================================================================
% MODEL PARAMETERS
% =========================================================================
ndim    = 3;            % number of spatial dimensions
nel     = MESH.nel;     % number of elements
EL2NOD  = uint32(MESH.EL2NOD{1});  % connectivity matrix for velocity problem
EL2NODP = uint32(MESH.EL2NOD{1}(1:4,:)); % connectivity matrix for pressure problem
nUdof   = ndim*max(max(EL2NOD));
nPdof   = max(max(EL2NODP));
nnodel  = size(EL2NOD,1);  % number of nodes in each element
nUdofel = ndim*nnodel;     % number of velocity dofs in each element
nPdofel = size(EL2NODP,1); % number of pressure dofs in each element


%==========================================================================
% PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
%==========================================================================
[IP_X,SF.IP_w]       = ip_tetrahedron(OPTS.nip);
    % local coordinates and weights of points for integration of
    % velocity/pressure matrices
[SF.NU,SF.dNUds]     = sf_dsf_tet(IP_X,nnodel,'cell');
    % velocity shape functions and their derivatives
[SF.NP,~]            = sf_dsf_tet(IP_X,nPdofel,'cell');
    % linear pressure shape functions and their derivatives
    % derivatives of linear (4-node) shape functions are used to calculate
    % each element's Jacobian
[~,dN4ds]            = sf_dsf_tet(IP_X,nPdofel,'cell');
SF.dN4ds             = dN4ds{1};
    % derivatives of linear (4-node) shape functions (used to calculate each element's Jacobian 
    % when the elements have straight edges). These are the same for any global coordinate system
[SF.NU14,SF.dNUds14] = sf_dsf_tet(IP_X,size(MESH.EL2NOD14,1),'cell');
    % shape functions and their derivatives for 14-node elements. They are 
    % used to calculate the Jacobian for isoparametric elements
[SF.NU20,SF.dNUds20] = sf_dsf_tet(IP_X,size(MESH.EL2NOD_cubic,1),'cell');
    % shape functions and their derivatives for 20-node elements. They are 
    % used to calculate the Jacobian for isoparametric elements

%==========================================================================
% STORAGE FOR GLOBAL MATRICES AND VECTORS
%==========================================================================
ASSEMBLY.KKv = zeros(nUdofel*(nUdofel+1)/2,nel); % storage for stiffness matrix data
ASSEMBLY.GGv = zeros(nPdofel*nUdofel      ,nel); % storage for gradient matrix data
ASSEMBLY.MMv = zeros(nPdofel*(nPdofel+1)/2,nel); % storage for (1/mu)-scaled mass matrix data
ASSEMBLY.Fb  = zeros(nUdofel,nel); % storage for global force vector
if isfield(VAR,'DilEl')
    ASSEMBLY.Fpdil = zeros(nPdofel,nel); % storage for global dilation vector (pressure eq.)
end

% Density and viscosity at integration points of all elements
if OPTS.return_dens_el
    ASSEMBLY.DensEl = zeros(nel,OPTS.nip);
end
if OPTS.return_visc_ip
    ASSEMBLY.ViscIP = zeros(nel,OPTS.nip);
end

%==========================================================================
% ASSEMBLY STIFFNESS MATRICES FOR EACH SET OF ELEMENTS
%==========================================================================
ASSEMBLY = assembly_stokes_els_out_cone_no_cross_2pi_opt...
    (ASSEMBLY,SF,VAR,MESH,PHYSICS,NUMSCALE,OPTS);
ASSEMBLY = assembly_stokes_els_out_cone_cross_2pi_opt...
    (ASSEMBLY,SF,VAR,MESH,PHYSICS,NUMSCALE,OPTS);
ASSEMBLY = assembly_stokes_els_in_cone_no_iso_opt...
    (ASSEMBLY,SF,VAR,MESH,PHYSICS,NUMSCALE,OPTS);
% ASSEMBLY = assembly_stokes_els_in_cone_iso_face_nodes_opt...
%     (ASSEMBLY,SF,VAR,MESH,PHYSICS,NUMSCALE,OPTS);
ASSEMBLY = assembly_stokes_els_in_cone_iso_20_nodes_opt...
            (ASSEMBLY,SF,VAR,MESH,PHYSICS,NUMSCALE,OPTS);

%==========================================================================
% ASSEMBLY OF GLOBAL STIFFNESS MATRIX
%==========================================================================
if OPTS.use_mutils
    opts_mutils.symmetric  = 1;
    opts_mutils.n_node_dof = ndim;
    if isunix
        opts_mutils.nthreads = OPTS.nthreads;
    else
        opts_mutils.nthreads = 1; % >1 is very slow (????)
    end
    KK = sparse_create(EL2NOD, ASSEMBLY.KKv, opts_mutils);
    
else
    % CREATE TRIPLET FORMAT INDICES
    EL2DOF               = zeros(nUdofel, nel,'int32');
    EL2DOF(1:ndim:end,:) = ndim*(EL2NOD-1)+1;
    EL2DOF(2:ndim:end,:) = ndim*(EL2NOD-1)+2;
    EL2DOF(3:ndim:end,:) = ndim*(EL2NOD-1)+3;
    indx_j     = repmat(int32(1:nUdofel),nUdofel,1); indx_i = indx_j';
    indx_i     = tril(indx_i); indx_i = indx_i(:); indx_i = indx_i(indx_i>0);
    indx_j     = tril(indx_j); indx_j = indx_j(:); indx_j = indx_j(indx_j>0);

    K_i        = EL2DOF(indx_i(:),:);
    K_j        = EL2DOF(indx_j(:),:);
    indx       = K_i < K_j;
    tmp        = K_j(indx);
    K_j(indx)  = K_i(indx);
    K_i(indx)  = tmp; clear indx tmp

    % convert triplet data to sparse global stiffness matrix (assembly)
    KK = sparse3(K_i , K_j , ASSEMBLY.KKv , nUdof , nUdof );
    clear K_i K_j
end
ASSEMBLY = rmfield(ASSEMBLY,'KKv');

%==========================================================================
% ASSEMBLY OF GLOBAL GRADIENT MATRIX
%==========================================================================
EL2DOF               = zeros(nUdofel,nel,'int32');
EL2DOF(1:ndim:end,:) = ndim*(EL2NOD-1)+1;
EL2DOF(2:ndim:end,:) = ndim*(EL2NOD-1)+2;
EL2DOF(3:ndim:end,:) = ndim*(EL2NOD-1)+3;
G_i    = repmat(EL2DOF,nPdofel,1); % global velocity dofs
indx_j = repmat(int32(1:nPdofel),nUdofel,1);
G_j    = EL2NODP(indx_j,:);        % global pressure nodes
% convert triplet data to sparse global gradient matrix (assembly)
GG     = sparse3(G_i , G_j , ASSEMBLY.GGv , nUdof , nPdof);
clear G_i G_j
ASSEMBLY = rmfield(ASSEMBLY,'GGv');

%==========================================================================
% ASSEMBLY OF GLOBAL FORCE VECTOR
%==========================================================================
Fb     = accumarray(EL2DOF(:) , ASSEMBLY.Fb(:));
clear EL2DOF
ASSEMBLY = rmfield(ASSEMBLY,'Fb');

%==========================================================================
% ASSEMBLY OF GLOBAL (1/VISCOSITY)-SCALED MASS MATRIX
%==========================================================================
if OPTS.use_mutils
    opts_mutils.symmetric  = 1;
    opts_mutils.n_node_dof = 1;
    if isunix
        opts_mutils.nthreads = OPTS.nthreads;
    else
        opts_mutils.nthreads = 1; % >1 is very slow (????)
    end
    MM = sparse_create(EL2NODP,ASSEMBLY.MMv,opts_mutils);
    
else
    indx_j    = repmat(int32(1:nPdofel),nPdofel,1); indx_i = indx_j';
    indx_i    = tril(indx_i); indx_i = indx_i(:); indx_i = indx_i(indx_i>0);
    indx_j    = tril(indx_j); indx_j = indx_j(:); indx_j = indx_j(indx_j>0);
    M_i       = EL2NODP(indx_i,:); M_i  = M_i(:);
    M_j       = EL2NODP(indx_j,:); M_j  = M_j(:);
    indx      = M_i < M_j;
    tmp       = M_j(indx);
    M_j(indx) = M_i(indx);
    M_i(indx) = tmp;
    % convert triplet data to sparse global matrix
    MM        = sparse3(M_i , M_j , ASSEMBLY.MMv , nPdof , nPdof);
end
clear M_i M_j
ASSEMBLY = rmfield(ASSEMBLY,'MMv');

%==========================================================================
% ASSEMBLY OF GLOBAL DILATION R.H.S VECTOR
%==========================================================================
if isfield(VAR,'DilEl')
    Fpdil = accumarray(EL2NODP(:),ASSEMBLY.Fpdil(:));
    clear EL2NODP
else
    Fpdil = 0;
end

% Density and viscosity at integration points of all elements
if OPTS.return_dens_el
    DensEl = ASSEMBLY.DensEl;
else
    DensEl = [];
end
if OPTS.return_visc_ip
    ViscIP = ASSEMBLY.ViscIP;
else
    ViscIP = [];
end
clear ASSEMBLY

end % END OF SUBFUNCTION element_assembly_opt