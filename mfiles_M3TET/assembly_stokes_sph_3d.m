function [KK,Fb,GG,MM,Fpdil,DensEl,ViscIP] = assembly_stokes_sph_3d...
    (VAR,MESH,PHYSICS,NUMSCALE,OPTS)
% Usage: [KK,Fb,GG,MM,Fpdil,DensEl,ViscIP] = assembly_stokes_sph_3d...
%   (VAR,MESH,PHYSICS,NUMSCALE,OPTS)
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
% JH Apr 2014 : cleaned up and made similar to visco-elastic assembly function
%

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

end % END OF FUNCTION assembly_stokes_sph_3d

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
GCOORD  = MESH.GCOORD;
PhaseID = MESH.PhaseID;
EL2NOD  = double(MESH.EL2NOD{1});
EL2NODP = double(MESH.EL2NOD{1}(1:4,:));
nnod    = MESH.nnod;
nel     = MESH.nel;
ndim    = 3;
nUdof   = nnod * ndim;
nPdof   = MESH.nVnod;
nnodel  = size(EL2NOD,1);
nUdofel = ndim*nnodel;
nPdofel = size(EL2NODP,1);


%==========================================================================
% PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
%==========================================================================
[IP_X,IP_w] = ip_tetrahedron(OPTS.nip);
    % local coordinates and weights of points for integration of
    % velocity/pressure matrices
[NU,dNUds]  = sf_dsf_tet(IP_X,nnodel,'cell');
    % velocity shape functions and their derivatives
[NP,dNPds]  = sf_dsf_tet(IP_X,nPdofel,'cell');
    % linear pressure shape functions and their derivatives
    % derivatives of linear (4-node) shape functions are used to calculate
    % each element's Jacobian
    
    
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
Fb    = zeros(nUdof,1);        % storage for global force vector
Fpdil = zeros(nPdof,1);        % storage for global dilation vector (pressure eq.)
KKv   = zeros(nUdofel*nUdofel,nel);  % storage for stiffness matrix data
KKi   = zeros(nUdofel*nUdofel,nel);  % storage for stiffness matrix i indices
KKj   = zeros(nUdofel*nUdofel,nel);  % storage for stiffness matrix j indices
GGv   = zeros(nPdofel*nUdofel,nel); % storage for gradient matrix data
GGi   = zeros(nPdofel*nUdofel,nel); % storage for gradient matrix i indices
GGj   = zeros(nPdofel*nUdofel,nel); % storage for gradient matrix j indices

% Inverse-(viscosity*penalty)-scaled Mass matrix for pressure preconditioning
MMv   = zeros(nPdofel*nPdofel,nel); % storage for (1/mu)-scaled mass matrix data
MMi   = zeros(nPdofel*nPdofel,nel); % storage for mass matrix i indices
MMj   = zeros(nPdofel*nPdofel,nel); % storage for mass matrix j indices

% Density and viscosity at integration points of all elements
if OPTS.return_dens_el
    DensEl = zeros(nel,OPTS.nip);
else
    DensEl = [];
end
if OPTS.return_visc_ip
    ViscIP = zeros(nel,OPTS.nip);
else
    ViscIP = [];
end


%==========================================================================
% ELEMENT TO GLOBAL VELOCITY D.O.F. CONNECTIVITY
%==========================================================================
EL2DOF                   = zeros(nUdofel, nel); % element -> global velocity dofs
EL2DOF(1:ndim:nUdofel,:) = ndim*EL2NOD-2;
EL2DOF(2:ndim:nUdofel,:) = ndim*EL2NOD-1;
EL2DOF(3:ndim:nUdofel,:) = ndim*EL2NOD;


%==========================================================================
% ELEMENT LOOP - MATRIX COMPUTATION
%==========================================================================
for iel = 1:nel
	% Initialize element arrays
    KKe = zeros(nUdofel);         % element stiffness matrix
    GGe = zeros(nUdofel,nPdofel); % element gradient matrix
    MMe = zeros(nPdofel);         % (1/mu)-scaled element mass matrix
    Fbe = zeros(nUdofel,1); % element force vector
    Fpe = zeros(nPdofel,1); % element dilatation vector (pressure eq.)

    if isfield(VAR,'DilEl')
        DilEl = VAR.DilEl(iel); % dilation prescribed for each element
    else
        DilEl = 0;
    end
    
    %======================================================================
    % INTEGRATION LOOP
    %======================================================================
    for ip=1:OPTS.nip
        %======================================================================
        % CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
        %======================================================================
        % NOTE: In spherical coordinates the elements may have curved
        %       faces/edges. Hence the Jacobian is different at every
        %       integration point and must be calculated inside the
        %       integration loop.
        ECOORD = GCOORD(:,EL2NOD(:,iel));
            % element node coordinates
        J      = dNUds{ip}'*ECOORD';
            % Jacobi matrix: reference element ==> current element
        detJ   =  J(1)*J(5)*J(9) + J(4)*J(8)*J(3) + J(7)*J(2)*J(6) ...
                - J(7)*J(5)*J(3) - J(1)*J(8)*J(6) - J(4)*J(2)*J(9);
            % Determinate of Jacobi matrix
        if(detJ<=0)
            err = ['negative jacobian in element ' iel];
            error(err);
        end
        
        xyz_ip  = ECOORD*NU{ip};          % physical coordinates of integration point
        r_ip    = sqrt(sum(xyz_ip.^2,1)); % radius of integration point
        rhat_ip = xyz_ip ./ r_ip;         % scale them all by l
        
        
        %==================================================================
        % PROPERTIES OF ELEMENTS AT ip-TH INTEGRATION POINT
        %==================================================================
        Dens_iel = calc_element_dens(VAR,OPTS,PHYSICS,EL2NOD,PhaseID,iel,NP{ip});
            % density, *SUBFUNCTION*
        Visc_iel = calc_element_visc(VAR,OPTS,PHYSICS,EL2NOD,PhaseID,iel,NP{ip});
            % viscosity; *SUBFUNCTION*
        
        if OPTS.return_dens_el
            DensEl(iel,ip) = Dens_iel(:);
        end
        if OPTS.return_visc_ip
            ViscIP(iel,ip) = Visc_iel(:);
        end
        
        % Gravitational force at ip-th integration point
        % Note that Bscale is new defined: Bscale = L0^2/(Visc0 * U0) 
        % It does no longer include grav. acceleration and density
        Fg_el = PHYSICS.g * NUMSCALE.Bscale * (Dens_iel-NUMSCALE.Dens0);

        weight = detJ * IP_w(ip);
            % element area * Gauss weight (weights includes volume of master
            % tetrahedron A=1/6) (V=0.5*A*h; A=area of base;h=height)

        dNUdx = J\dNUds{ip}'; % global derivative of shape function

        B(1,1:3:nUdofel-2) = dNUdx(1,:); % Hughes' book p.87 2.8.21 (here transpose of B)
        B(2,2:3:nUdofel-1) = dNUdx(2,:); % 2.8.21 shows the block for one node
        B(3,3:3:nUdofel  ) = dNUdx(3,:); % here: 10 blocks after another
        B(4,2:3:nUdofel-1) = dNUdx(3,:); B(4,3:3:nUdofel  ) = dNUdx(2,:);
        B(5,1:3:nUdofel-2) = dNUdx(3,:); B(5,3:3:nUdofel  ) = dNUdx(1,:);
        B(6,1:3:nUdofel-2) = dNUdx(2,:); B(6,2:3:nUdofel-1) = dNUdx(1,:);
        
        De  = D * Visc_iel;

        % Assemble element matrices
        KKe = KKe  + (B'*De*B) * weight;
            % element stiffness matrix
        GGe = GGe  + (dNUdx(:)*NP{ip}') * weight;
            % element gradient matrix
        MMe = MMe  + (NP{ip}*NP{ip}') * (weight./Visc_iel);
            % (1/viscosity)-scaled element mass matrix
        
        % Assemble force vector
%         Fbe(3:3:nUdofel) = Fbe(3:3:nUdofel) + NU{ip} .* Fg_el * weight; % buoyancy force for cartesian coords
        Fbe(1:3:nUdofel) = Fbe(1:3:nUdofel) + NU{ip} .* Fg_el * rhat_ip(1) * weight; % buoyancy force x-dir
        Fbe(2:3:nUdofel) = Fbe(2:3:nUdofel) + NU{ip} .* Fg_el * rhat_ip(2) * weight; % buoyancy force y-dir
        Fbe(3:3:nUdofel) = Fbe(3:3:nUdofel) + NU{ip} .* Fg_el * rhat_ip(3) * weight; % buoyancy force z-dir

        % R.h.s. term on incompressibility equation from non-zero dilation
        Fpe              = Fpe + NP{ip}*DilEl*weight;
    end % END OF INTEGRATION LOOP
    
    %======================================================================
    % STORE DATA OF ALL ELEMENTS IN BLOCK FOR ASSEMBLY
    %======================================================================
    gUdof = EL2DOF(:,iel);   % global velocity dofs of current element
    gPdof = EL2NOD(1:4,iel); % global pressure dofs of current element
    
    % Store stiffness matrix of each element for assembly
    rows        = gUdof*ones(1,nUdofel);
    cols        = ones(nUdofel,1)*gUdof';
    KKv(:,iel)  = KKe(:);  % values
    KKi(:,iel)  = rows(:); % i indices
    KKj(:,iel)  = cols(:); % j indices

    % Assemble global force vector including body force term
    Fb(gUdof)   = Fb(gUdof) + Fbe;
    
    % Assemble global pressure-dilation r.h.s. vector
    Fpdil(gPdof)= Fpdil(gPdof) + Fpe;
    
    % Store gradient matrix of each element for assembly
    rows        = gUdof*ones(1,nPdofel);  % rows of GGe are velocity dofs
    cols        = ones(nUdofel,1)*gPdof'; % cols of GGe are pressure dofs 
    GGv(:,iel)  = GGe(:);  % values
    GGi(:,iel)  = rows(:); % i indices
    GGj(:,iel)  = cols(:); % j indices
    
    % Store (1/viscosity)-scaled mass matrix of each element
    rows        = gPdof*ones(1,nPdofel);
    cols        = ones(nPdofel,1)*gPdof';
    MMv(:,iel)  = MMe(:);  % values
    MMi(:,iel)  = rows(:); % i indices
    MMj(:,iel)  = cols(:); % j indices
end % END OF ELEMENT LOOP

% Average density in each element is sufficient
if OPTS.return_dens_el
    DensEl = mean(DensEl,2);
end

%==========================================================================
% ASSEMBLY OF GLOBAL SPARSE MATRICES AND RHS-VECTOR
%==========================================================================
MM = sparse3(MMi(:),MMj(:),MMv(:)); % global (1/viscosity)-scaled mass matrix
GG = sparse3(GGi(:),GGj(:),GGv(:)); % global gradient matrix
KK = sparse3(KKi(:),KKj(:),KKv(:)); % global stiffness matrix
clear GGi GGj GGv MMi MMj MMv KKi KKj KKv

% only keep lower half for compatibility with optimized version
MM = tril(MM);
KK = tril(KK);

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
GCOORD  = MESH.GCOORD; % node coordinates
PhaseID = MESH.PhaseID;
if length(PhaseID)==1; PhaseID = ones(MESH.nel,1); end
nUdof   = ndim*max(max(EL2NOD));
nPdof   = max(max(EL2NODP));
nnodel  = size(EL2NOD,1);  % number of nodes in each element
nvertx  = MESH.nvertx;     % number of vertices (or corners) in an element
nUdofel = ndim*nnodel;     % number of velocity dofs in each element
nPdofel = size(EL2NODP,1); % number of pressure dofs in each element


%==========================================================================
% PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
%==========================================================================
[IP_X,IP_w] = ip_tetrahedron(OPTS.nip);
    % local coordinates and weights of points for integration of
    % velocity/pressure matrices
[NU,dNUds]  = sf_dsf_tet(IP_X,nnodel,'cell');
    % velocity shape functions and their derivatives
[NP,dNPds]  = sf_dsf_tet(IP_X,nPdofel,'cell');
    % linear pressure shape functions and their derivatives
    % derivatives of linear (4-node) shape functions are used to calculate
    % each element's Jacobian


%==========================================================================
% STORAGE FOR DATA OF ALL ELEMENT MATRICES/VECTORS
%==========================================================================
Fb_all  = zeros(nUdofel              , nel);  % storage for global force vector
Fp_all  = zeros(nPdofel              , nel);  % storage for global dilation force vector
K_all   = zeros(nUdofel*(nUdofel+1)/2, nel);  % storage for stiffness matrix data
G_all   = zeros(nPdofel*nUdofel      , nel);  % storage for gradient matrix data
M_all   = zeros(nPdofel*(nPdofel+1)/2, nel);  % storage for (1/mu)-scaled global mass matrix

% Density and viscosity at integration points of all elements
if OPTS.return_dens_el
    DensEl = zeros(nel,OPTS.nip);
else
    DensEl = [];
end
if OPTS.return_visc_ip
    ViscIP = zeros(nel,OPTS.nip);
else
    ViscIP = [];
end


%==========================================================================
% BLOCKING PARAMETERS (nelblk must be < nel)
%==========================================================================
nelblk    = min(nel,OPTS.nelblk); % in case nel<nelblk
nblk      = ceil(nel/nelblk);     % number of blocks
nel_asmbl = min(nel,OPTS.nblk_asmbl*nelblk); % number of elements per partial K assembly
il        = 1;
iu        = nelblk;
ilK       = 1;      % first element in K-block
iuK       = nelblk; % last element in K-block
el1K      = 1;                  % global number of first element in K-block
el2K      = min(nel_asmbl,nel); % global number of last element in K-block

tic
%==========================================================================
% BLOCK LOOP - MATRIX COMPUTATION
%==========================================================================
for ib=1:nblk % loop over element blocks
    %======================================================================
    % STORAGE FOR DATA OF ELEMENTS IN BLOCK
    %======================================================================
    K_blk  = zeros(nelblk, nUdofel*(nUdofel+1)/2);
        % symmetric stiffness matrix, dim=vel dofs, but only upper triangle
    M_blk  = zeros(nelblk, nPdofel*(nPdofel+1)/2);
        % symmetric (1/mu)-scaled mass matrix, dim=pressure nodes, but only upper triangle
    G_blk  = zeros(nelblk, nPdofel*nUdofel);
        % asymmetric gradient matrix, vel dofs x pressure nodes
    Fb_blk = zeros(nelblk, nUdofel);
        % storage for entries in bouyancy force vector
    Fp_blk = zeros(nelblk, nPdofel);
        % storage for entries in dilation rhs vector

    if isfield(VAR,'DilEl')
        Dil_blk = VAR.DilEl(il:iu); % dilation prescribed for each element
    else
        Dil_blk = [];
    end
    
    %======================================================================
    % INTEGRATION LOOP: DEAL WITH A BLOCK OF ELEMENTS AT ONCE
    %======================================================================
    for ip=1:OPTS.nip
        %======================================================================
        % CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
        %======================================================================
        % NOTE: In spherical coordinates the elements may have curved
        %       faces/edges. Hence the Jacobian is different at every
        %       integration point and must be calculated inside the
        %       integration loop.
        
        % COORDINATES OF CORNER/VERTEX NODES OF LL ELEMENTS IN BLOCK
        ECOORD_x = reshape( GCOORD(1,EL2NOD(:,il:iu)), nnodel, nelblk );
        ECOORD_y = reshape( GCOORD(2,EL2NOD(:,il:iu)), nnodel, nelblk );
        ECOORD_z = reshape( GCOORD(3,EL2NOD(:,il:iu)), nnodel, nelblk );
        Jx       = ECOORD_x'*dNUds{ip};
        Jy       = ECOORD_y'*dNUds{ip};
        Jz       = ECOORD_z'*dNUds{ip};

        detJ     = elblk_detA(Jx,Jy,Jz);
        
        [invJx,invJy,invJz] = elblk_invert_A(Jx,Jy,Jz,detJ);
        
        %==================================================================
        % DERIVATIVES wrt GLOBAL COORDINATES
        %==================================================================
        dNUdx    = invJx*dNUds{ip}';
        dNUdy    = invJy*dNUds{ip}';
        dNUdz    = invJz*dNUds{ip}';
        
        % physical coordinates of integration point
        x_ip_blk = ECOORD_x' * NU{ip};
        y_ip_blk = ECOORD_y' * NU{ip};
        z_ip_blk = ECOORD_z' * NU{ip};
        
        % radius of integration point
        r_ip_blk = sqrt( x_ip_blk.^2 + y_ip_blk.^2 + z_ip_blk.^2 );
        
        % scale them all by l
        xhat_ip_blk = x_ip_blk ./ r_ip_blk;
        yhat_ip_blk = y_ip_blk ./ r_ip_blk;
        zhat_ip_blk = z_ip_blk ./ r_ip_blk;
        
        % Integration points in spherical coordinates
        [ip_coord_sph] = cartesian2spherical([x_ip_blk'; y_ip_blk'; z_ip_blk']);
        th_ip          = ip_coord_sph(1,:)';
        ph_ip          = ip_coord_sph(2,:)';
        r_ip           = ip_coord_sph(3,:)';
        if max(abs(r_ip-r_ip_blk))>1e-12
            error('Transformation of cartesian to spherical coordinates failed.');
        end
        
        %==================================================================
        % PROPERTIES OF ELEMENTS AT ip-TH INTEGRATION POINT
        %==================================================================
        Dens_blk = calc_element_dens(VAR,OPTS,PHYSICS,EL2NOD,PhaseID,il:iu,NP{ip},th_ip,ph_ip,r_ip);
            % density, *SUBFUNCTION*
        Visc_blk = calc_element_visc(VAR,OPTS,PHYSICS,EL2NOD,PhaseID,il:iu,NP{ip},th_ip,ph_ip,r_ip);
            % viscosity; *SUBFUNCTION*
        if OPTS.return_dens_el
            DensEl(il:iu,ip) = Dens_blk(:);
        end
        if OPTS.return_visc_ip
            ViscIP(il:iu,ip) = Visc_blk(:);
        end
        
        % Gravitational force at ip-th integration point
        % Note that Bscale is new defined: Bscale = L0^2/(Visc0 * U0) 
        % (it does no longer includes grav. acceleration and density)
        Fg_blk = PHYSICS.g * NUMSCALE.Bscale .* (Dens_blk-NUMSCALE.Dens0);

        %==================================================================
        % NUMERICAL INTEGRATION OF ELEMENT MATRICES
        %==================================================================
        weight  = detJ*IP_w(ip); % integration weight times volumes of tetrahedrons
        
        % ASSEMBLY OF STIFFNESS MATRICES FOR ALL ELEMENTS IN BLOCK
        switch OPTS.type_Dmat
            case '23rd'
                C1 =  4/3; % Used instead of 'D' matrix in standard assembly
                C2 = -2/3; % see Zienkiewicz book, Vol 2, 4th edition, p 519
            case '221'
                C1 = 2;
                C2 = 0;
        end
        weight2 = Visc_blk.*weight;
            % "weight" times viscosity at integration point in all elements
        indx    = 1;
        if C2==0 % AVOID MULTIPLICATIONS BY ZERO
            for i=1:nnodel
                % x-velocity equation (1st, 4th, 7th,... rows of stiffness matrices)
                for j=i:nnodel
                    % x-velocity (1st, 4th, 7th,... columns)
                    K_blk(:,indx) = K_blk(:,indx) + weight2 .* ...
                        (C1*dNUdx(:,i).*dNUdx(:,j) + dNUdy(:,i).*dNUdy(:,j) + dNUdz(:,i).*dNUdz(:,j));
                    indx = indx + 1;
                    % y-velocity (2nd, 5th, 8th,... columns)
                    K_blk(:,indx) = K_blk(:,indx) + weight2 .* (dNUdy(:,i).*dNUdx(:,j));
                    indx = indx + 1;
                    % z-velocity (3rd, 6th, 9th,... columns)
                    K_blk(:,indx) = K_blk(:,indx) + weight2 .* (dNUdz(:,i).*dNUdx(:,j));
                    indx = indx + 1;
                end

                % y-velocity equation (2nd, 5th, 8th,... rows of stiffness matrices)
                for j=i:nnodel
                    if j>i
                        %  x-velocity (4th, 7th, 10th,... columns)
                        K_blk(:,indx) = K_blk(:,indx) + weight2 .* (dNUdx(:,i).*dNUdy(:,j));
                        indx = indx + 1;
                    end
                    % y-velocity (2nd, 5th, 8th,... columns)
                    K_blk(:,indx) = K_blk(:,indx) + weight2 .* ...
                        (C1*dNUdy(:,i).*dNUdy(:,j) + dNUdx(:,i).*dNUdx(:,j) + dNUdz(:,i).*dNUdz(:,j));
                    indx = indx + 1;
                    % z-velocity (3rd, 6th, 9th,... columns)
                    K_blk(:,indx) = K_blk(:,indx) + weight2 .* (dNUdz(:,i).*dNUdy(:,j));
                    indx = indx + 1;
                end

                % z-velocity equation (3rd, 6th, 9th,... rows of stiffness matrices)
                for j=i:nnodel
                    if j>i
                        %  x-velocity (4th, 7th, 10th,... columns)
                        K_blk(:,indx) = K_blk(:,indx) + weight2 .* (dNUdx(:,i).*dNUdz(:,j));
                        indx = indx + 1;

                        % y-velocity (5th, 8th, 11th,... columns)
                        K_blk(:,indx) = K_blk(:,indx) + weight2 .* (dNUdy(:,i).*dNUdz(:,j));
                        indx = indx + 1;
                    end
                    % z-velocity (3rd, 6th, 9th,... columns)
                    K_blk(:,indx) = K_blk(:,indx) + weight2 .* ...
                        (C1*dNUdz(:,i).*dNUdz(:,j) + dNUdx(:,i).*dNUdx(:,j) + dNUdy(:,i).*dNUdy(:,j));
                    indx = indx + 1;
                end
            end
            
        else
            for i=1:nnodel
                % x-velocity equation (1st, 4th, 7th,... rows of stiffness matrices)
                for j=i:nnodel
                    % x-velocity (1st, 4th, 7th,... columns)
                    K_blk(:,indx) = K_blk(:,indx) + weight2 .* ...
                        (C1*dNUdx(:,i).*dNUdx(:,j) + dNUdy(:,i).*dNUdy(:,j) + dNUdz(:,i).*dNUdz(:,j));
                    indx = indx + 1;
                    % y-velocity (2nd, 5th, 8th,... columns)
                    K_blk(:,indx) = K_blk(:,indx) + weight2 .* ...
                        (C2*dNUdx(:,i).*dNUdy(:,j) + dNUdy(:,i).*dNUdx(:,j));
                    indx = indx + 1;
                    % z-velocity (3rd, 6th, 9th,... columns)
                    K_blk(:,indx) = K_blk(:,indx) + weight2 .* ...
                        (C2*dNUdx(:,i).*dNUdz(:,j) + dNUdz(:,i).*dNUdx(:,j));
                    indx = indx + 1;
                end

                % y-velocity equation (2nd, 5th, 8th,... rows of stiffness matrices)
                for j=i:nnodel
                    if j>i
                        %  x-velocity (4th, 7th, 10th,... columns)
                        K_blk(:,indx) = K_blk(:,indx) + weight2 .* ...
                            (C2*dNUdy(:,i).*dNUdx(:,j) + dNUdx(:,i).*dNUdy(:,j));
                        indx = indx + 1;
                    end
                    % y-velocity (2nd, 5th, 8th,... columns)
                    K_blk(:,indx) = K_blk(:,indx) + weight2 .* ...
                        (C1*dNUdy(:,i).*dNUdy(:,j) + dNUdx(:,i).*dNUdx(:,j) + dNUdz(:,i).*dNUdz(:,j));
                    indx = indx + 1;
                    % z-velocity (3rd, 6th, 9th,... columns)
                    K_blk(:,indx) = K_blk(:,indx) + weight2 .* ...
                        (C2*dNUdy(:,i).*dNUdz(:,j) + dNUdz(:,i).*dNUdy(:,j));
                    indx = indx + 1;
                end

                % z-velocity equation (3rd, 6th, 9th,... rows of stiffness matrices)
                for j=i:nnodel
                    if j>i
                        %  x-velocity (4th, 7th, 10th,... columns)
                        K_blk(:,indx) = K_blk(:,indx) + weight2 .* ...
                            (C2*dNUdz(:,i).*dNUdx(:,j) + dNUdx(:,i).*dNUdz(:,j));
                        indx = indx + 1;

                        % y-velocity (5th, 8th, 11th,... columns)
                        K_blk(:,indx) = K_blk(:,indx) + weight2 .* ...
                            (C2*dNUdz(:,i).*dNUdy(:,j) + dNUdy(:,i).*dNUdz(:,j));
                        indx = indx + 1;
                    end
                    % z-velocity (3rd, 6th, 9th,... columns)
                    K_blk(:,indx) = K_blk(:,indx) + weight2 .* ...
                        (C1*dNUdz(:,i).*dNUdz(:,j) + dNUdx(:,i).*dNUdx(:,j) + dNUdy(:,i).*dNUdy(:,j));
                    indx = indx + 1;
                end
            end
        end
        
        % ASSEMBLY OF GRADIENT MATRICES FOR ALL ELEMENTS IN BLOCK
        NP_blk = repmat(NP{ip}',nelblk,1);
        for i=1:nPdofel
            tmp1 = weight.*NP_blk(:,i);
            tmp2 = tmp1(:,ones(1,nnodel));
            ii   = (i-1)*nUdofel + (1:3:nUdofel-2);
            G_blk(:,ii) = G_blk(:,ii) + tmp2.*dNUdx;
            ii   = (i-1)*nUdofel + (2:3:nUdofel-1);
            G_blk(:,ii) = G_blk(:,ii) + tmp2.*dNUdy;
            ii   = (i-1)*nUdofel + (3:3:nUdofel);
            G_blk(:,ii) = G_blk(:,ii) + tmp2.*dNUdz;
        end

        % ASSEMBLY OF (1/VISCOSITY)-SCALED MASS MATRICES FOR ALL ELEMENTS IN BLOCK
        weight3 = weight./Visc_blk;
            % "weight" divided by viscosity at integration point in all elements
        indx    = 1;
        for i=1:nPdofel
            for j=i:nPdofel
                M_blk(:,indx) = M_blk(:,indx) + weight3 .* ...
                                NP{ip}(i)*NP{ip}(j);
                                % ==NP_blk(:,i) .* NP_blk(:,j);
                indx = indx + 1;
            end
        end
        
        % ASSEMBLY OF BUOYANCY FORCE VECTORS FOR ALL ELEMENTS IN BLOCK
        Fb_blk(:,1:3:nUdofel) = Fb_blk(:,1:3:nUdofel) + (xhat_ip_blk .* weight .* Fg_blk) * NU{ip}';
        Fb_blk(:,2:3:nUdofel) = Fb_blk(:,2:3:nUdofel) + (yhat_ip_blk .* weight .* Fg_blk) * NU{ip}';
        Fb_blk(:,3:3:nUdofel) = Fb_blk(:,3:3:nUdofel) + (zhat_ip_blk .* weight .* Fg_blk) * NU{ip}';
        
        % ASSEMBY OF DILATATION RHS VECTORS FOR ALL ELEMENTS IN BLOCK
        % (a 'force'-like term on the incompressibility equation results
        % from a non-zero dilation)
        if ~isempty(Dil_blk)
            Fp_blk = Fp_blk + (const .* Dil_blk) * NP{ip}';
        end
    end % END OF INTEGRATION LOOP
    
    %======================================================================
    % STORE DATA OF ALL ELEMENTS IN BLOCK FOR ASSEMBLY
    %======================================================================
    K_all(:,ilK:iuK) = K_blk';
    G_all(:,il:iu)   = G_blk';
    M_all(:,il:iu)   = M_blk';
    Fb_all(:,il:iu)  = Fb_blk';
    if ~isempty(Dil_blk)
        Fp_all(:,il:iu)  = Fp_blk';
    end
    
    %======================================================================
    % PARTIAL ASSEMBLY OF GLOBAL MATRICES TO AVOID PEAKS IN MEMORY USAGE
    %======================================================================
    if OPTS.use_mutils && ib==nblk
        opts_mutils.symmetric  = 1;
        opts_mutils.n_node_dof = ndim;
        if isunix
            opts_mutils.nthreads = OPTS.nthreads;
        else
            opts_mutils.nthreads = 1; % >1 is very slow (????)
        end
        K_all(abs(K_all)<OPTS.rtol_asmbly*max(abs(K_all(:)))) = 0;
        KK = sparse_create(EL2NOD, K_all, opts_mutils);
        clear K_all
    end
    if ~OPTS.use_mutils && (iuK==nel_asmbl || ib==nblk)
        if ib==nblk
            K_all(:,iuK+1:end) = [];
            nel_asmbl          = el2K-el1K+1;
        end
        
        %==================================================================
        % CREATE TRIPLET FORMAT INDICES
        %==================================================================
        EL2DOF               = zeros(nUdofel, nel_asmbl,'int32');
        EL2DOF(1:ndim:end,:) = ndim*(EL2NOD(:,el1K:el2K)-1)+1;
        EL2DOF(2:ndim:end,:) = ndim*(EL2NOD(:,el1K:el2K)-1)+2;
        EL2DOF(3:ndim:end,:) = ndim*(EL2NOD(:,el1K:el2K)-1)+3;
        indx_j     = repmat(int32(1:nUdofel),nUdofel,1); indx_i = indx_j';
        indx_i     = tril(indx_i); indx_i = indx_i(:); indx_i = indx_i(indx_i>0);
        indx_j     = tril(indx_j); indx_j = indx_j(:); indx_j = indx_j(indx_j>0);
        
        %==================================================================
        % PARTIAL ASSEMBLY OF GLOBAL STIFFNESS MATRIX
        %==================================================================
        K_i        = EL2DOF(indx_i(:),:);
        K_j        = EL2DOF(indx_j(:),:);
        indx       = K_i < K_j;
        tmp        = K_j(indx);
        K_j(indx)  = K_i(indx);
        K_i(indx)  = tmp; clear indx tmp
%         byte2MB = 1/1048576;
%         var = whos('K'); K_MB = var.bytes * byte2MB;
%         var = whos('K_i'); Ki_MB = var.bytes * byte2MB;
%         var = whos('K_j'); Kj_MB = var.bytes * byte2MB;
%         var = whos('K_all'); Kall_MB = var.bytes * byte2MB;

        % convert triplet data to sparse global stiffness matrix (assembly)
        ind = abs(K_all)>OPTS.rtol_asmbly*max(abs(K_all(:)));
        if el1K==1 % first partial assembly
            KK = sparse3(K_i(ind) , K_j(ind) , K_all(ind) , nUdof , nUdof );
        else
            KK = KK + sparse3(K_i(ind) , K_j(ind) , K_all(ind) , nUdof , nUdof );
        end
        
        if el2K==nel
            % Done with K-matrix assembly
            clear indx_i indx_j K_i K_j K_all
            break
        end
        el1K = el2K+1;
        el2K = min(nel,el2K+nel_asmbl);
        ilK  = 1; % reset counters for elements in K-block assembly
        iuK  = 0; %
    else
        ilK = ilK + nelblk; % update counter for elements in K-block assembly
    end
    
    %======================================================================
    % READJUST START, END AND SIZE OF BLOCK
    %======================================================================
    il = il + nelblk;
    if ib==nblk-1
        nelblk   = nel-iu;
    end
    iu  = iu  + nelblk; % update counter for elements in block
    iuK = iuK + nelblk; % update counter for elements in K-block
    
end % END OF BLOCK LOOP

% Done with accumulating element blocks. All data now in M_all,
% G_all, etc
clear K_blk M_blk G_blk Fb_blk Fp_blk invJx invJy invJz
            
% Average density in each element is sufficient
if OPTS.return_dens_el
    DensEl = mean(DensEl,2);
end

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
ind    = abs(G_all)>OPTS.rtol_asmbly*max(abs(G_all(:)));
GG     = sparse3(G_i(ind) , G_j(ind) , G_all(ind) , nUdof , nPdof);
clear G_i G_j G_all

%==========================================================================
% ASSEMBLY OF GLOBAL FORCE VECTOR
%==========================================================================
Fb     = accumarray(EL2DOF(:) , Fb_all(:));
clear Fb_all EL2DOF

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
    M_all(abs(M_all)<OPTS.rtol_asmbly*max(abs(M_all(:)))) = 0;
    MM = sparse_create(EL2NODP,M_all,opts_mutils);
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
    ind       = abs(M_all)>OPTS.rtol_asmbly*max(abs(M_all(:)));
    MM        = sparse3(M_i(ind) , M_j(ind) , M_all(ind) , nPdof , nPdof);
end
clear M_i M_j M_all

%==========================================================================
% ASSEMBLY OF GLOBAL DILATION R.H.S VECTOR
%==========================================================================
if isempty(Dil_blk)
    Fpdil = 0;
else
    Fpdil   = accumarray(EL2NODP(:),Fp_all(:));
    clear Fp_all EL2NODP
end

end % END OF FUNCTION element_assembly_opt