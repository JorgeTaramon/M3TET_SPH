function Vars_nod = ipval_to_nodval_3d(GCOORD,EL2NOD,Vars_ip,BC)
% Usage: Vars_nod = ipval_to_nodval_3d(GCOORD,EL2NOD,Vars_ip,BC)
% 
% Purpose: Calculates nodal values of one (or more) variables evaluated at
%          integration points. Makes a continuous, least-squares solution.
%
% Input:
%   GCOORD : [matrix]    : coordinates of all nodes in mesh
%   EL2NOD : [matrix]    : finite element connectivity matrix (nnodel x nel)
%   Vars_ip: [cell]      : variable values at integration points; format:
%            a cell of length "nvar" containing matrices of size nel x nip
%   BC     : [structure] : Boundary conditions on the variable (optional)
%
% Output:
%   Vars_nod : [cell]    : continuous variable fields at nodes; format:
%            a cell of length "nvar" containing vectors of size nnod x 1
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Oct 2016
%


method_assembly = 'opt';  % std=standard assembly (slow); opt=optimized (fast)
nelblk          = 1000;   % blocksize (for method_assembly='opt')
method_solve_eq = 'chol'; % 'chol' - Cholesky factorization
                          % 'pcg'  - preconditioned conjugate gradient
                          % '\'    - Matlab's backslash
return_format = 'cell';
if ~iscell(Vars_ip)
    nel = size(EL2NOD,2);
    if size(Vars_ip)~=nel
        error('Variable(s) must be provided as a cell or a single matrix of size (nel x nip)');
    end
    return_format = 'vector';
    Vars_ip       = {Vars_ip};
end

% Numebr of variables at integration points
nvar = length(Vars_ip);

% number of integration points
nip  = size(Vars_ip{1},2);

% local coordinates and weights of integration points
[IP_X,IP_w] = ip_tetrahedron(nip);

% shape functions and their derivatives        
[N,dNds]    = sf_dsf_tet(IP_X,4,'cell');

if nargin>3 && ~isempty(BC) && isfield(BC,'nfix') && BC.nfix>0
    no_BC = 0;
else
    no_BC = 1;
end

nnod     = max(max(EL2NOD(1:4,:)));
Vars_nod = cell(1,nvar);

% Assemble mass matrix MM
% =======================
switch method_assembly
    case 'std' % standard element assembly: good for understanding FE
        MM = assembly_mass_matrix_std(GCOORD,EL2NOD(1:4,:),N,dNds,IP_w); % *SUBFUNCTION*

    case 'opt' % a faster block-wise assembly (MILAMIN style)
        MM = assembly_mass_matrix_opt(GCOORD,EL2NOD(1:4,:),N,dNds,IP_w,nelblk); % *SUBFUNCTION*
end

for ivar=1:nvar
    % Integrate variable (calculate rhs of equation)
    % ==============================================
    switch method_assembly
        case 'std' % standard element assembly: good for understanding FE
            rhs = assembly_rhs_std(GCOORD,EL2NOD(1:4,:),N,dNds,IP_w,Vars_ip{ivar}); % *SUBFUNCTION*

        case 'opt' % a faster block-wise assembly (MILAMIN style)
            nelblk = 1000;
            rhs = assembly_rhs_opt(GCOORD,EL2NOD(1:4,:),N,dNds,IP_w,Vars_ip{ivar},nelblk); % *SUBFUNCTION*
    end

    % Solve equation to obtain nodal values
    % =====================================
    if no_BC
        switch method_solve_eq
            case 'pcg'
                % X = PCG(A,B,TOL,MAXIT,M1,M2,X0)
                V_nod = pcg(MM,rhs,1e-8,1000);
            otherwise % Cannot use Cholesky factorization in this case
                V_nod = MM \ rhs;
        end
        
    else
        ifix        = BC.ifix;
        vfix        = BC.vfix;
        V_nod       = zeros(nnod,1);
        V_nod(ifix) = vfix; % write prescribed values
        rhs         = rhs - MM(:,ifix) * V_nod(ifix,ivar); % move BC to rhs
        ifree       = setdiff(1:nnod,vfix);
        rhs         = rhs(ifree);
        
        switch method_solve_eq
            case 'chol' % CHOLESKY DIRECT SOLVER
                if ivar==1
                    % CHOLESKY FACTORIZATION
                    [L,perm] = cholesky_factorization(MM(ifree,ifree));
                end
                % FORWARD & BACKWARD SOLVE of triangular system
                V_nod(ifree) = cholesky_solve(rhs,L,perm);
                
            case 'pcg' % CONJUGATE GRADIENT SOLVER
                V_nod(ifree) = pcg(MM(ifree,ifree),rhs,1e-8,1000);
                
            case '\' % MATLAB'S BACKSLASH 
                V_nod(ifree) = MM(ifree,ifree) \ rhs;
        end
        
    end
    
    if size(EL2NOD,1)==10
        tmp           = V_nod;
        V_nod         = zeros(max(EL2NOD(:)),1);
        V_nod(1:nnod) = tmp;
        V_nod(EL2NOD(5:10,:)) = 0.5*(tmp(EL2NOD([1 2 3 4 1 2],:)) + ...
                                     tmp(EL2NOD([2 3 4 1 3 4],:)));
        clear tmp
    end
    
    if strcmp(return_format,'vector')
        Vars_nod = V_nod;
    else
        Vars_nod{ivar} = V_nod;
    end
end

end % END OF FUNCTION ipval_to_nodval_3d

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function MM = assembly_mass_matrix_std(GCOORD,EL2NOD,N,dNds,IP_w)

nnodel = size(EL2NOD,1);
nel    = size(EL2NOD,2);
Mv     = zeros(nnodel*nnodel,nel); % storage for mass matrix coefficients
iMC    = zeros(nnodel*nnodel,nel); % storage for the matrix i indixes
jMC    = zeros(nnodel*nnodel,nel); % storage for the matrix j indixes
nip    = length(IP_w);

for iel = 1:nel % main element loop
    % initialize element arrays
    elnod  = double(EL2NOD(1:4,iel))'; % element nodes
    ecoord = GCOORD(:,elnod); % coordinates of vertex nodes
    MMel   = zeros(nnodel);    % element capacity matrix
    for ip=1:nip % Integration loop
        % Jacobi matrix (mapping reference element ==> current element)
        J      = dNds{ip}'*ecoord'; 
        % Determinant of Jacobi matrix
        detJ   =  J(1)*J(5)*J(9) + J(4)*J(8)*J(3) + J(7)*J(2)*J(6) ...
                - J(7)*J(5)*J(3) - J(1)*J(8)*J(6) - J(4)*J(2)*J(9);
%         detJ  = det(J); % checking
        weight = detJ*IP_w(ip);
        
        % calculate element matrices
        MMel  = MMel + (N{ip}*N{ip}') * weight;
    end % integration
    
    % Store conductivity and capacity matrix of each element for assembly
    rows       = elnod'*ones(1,nnodel); % row location in global matrices
    cols       = ones(nnodel,1)*elnod;  % column location in global matrices
    Mv(:,iel)  = MMel(:);
    iMC(:,iel) = rows(:);
    jMC(:,iel) = cols(:);
end % element loop

% Assemble sparse matrix
MM  = sparse3(iMC(:),jMC(:),Mv(:));

end % END OF SUBFUNCTION assembly_mass_matrix_std

% #########################################################################

function NV = assembly_rhs_std(GCOORD,EL2NOD,N,dNds,IP_w,V_ip)

nnodel = size(EL2NOD,1);
nel    = size(EL2NOD,2);
NVv    = zeros(nnodel,nel);
nip    = length(IP_w);

for iel = 1:nel % main element loop
    % initialize element arrays
    elnod  = double(EL2NOD(1:4,iel))'; % element nodes
    ecoord = GCOORD(:,elnod); % coordinates of vertex nodes
    NVel   = zeros(nnodel,1);
    for ip=1:nip % Integration loop
        % Jacobi matrix (mapping reference element ==> current element)
        J      = dNds{ip}'*ecoord'; 
        % Determinant of Jacobi matrix
        detJ   =  J(1)*J(5)*J(9) + J(4)*J(8)*J(3) + J(7)*J(2)*J(6) ...
                - J(7)*J(5)*J(3) - J(1)*J(8)*J(6) - J(4)*J(2)*J(9);
%         detJ  = det(J); % checking
        weight = detJ*IP_w(ip);
        NVel   = NVel + N{ip}*V_ip(iel,ip) * weight;
    end % integration
    NVv(:,iel) = NVel;
end % element loop

% Assemble vector
iNV = EL2NOD';
NV  = accumarray(iNV(:),NVv(:));

end % END OF SUBFUNCTION assembly_rhs_std

% #########################################################################

function MM = assembly_mass_matrix_opt(GCOORD,EL2NOD,N,dNds,IP_w,nelblk)

nnodel = size(EL2NOD,1);
nel    = size(EL2NOD,2);
nelblk = min(nel,nelblk);
nblk   = ceil(nel/nelblk);
nip    = length(IP_w);

M_blk  = zeros(nelblk,nnodel*(nnodel+1)/2);
M_all  = zeros(nnodel*(nnodel+1)/2,nel);

il = 1;
iu = nelblk;
for iblk=1:nblk % loop over element blocks
    x_el     = reshape( GCOORD(1,EL2NOD(1:4,il:iu,:)), nnodel, nelblk );
    y_el     = reshape( GCOORD(2,EL2NOD(1:4,il:iu,:)), nnodel, nelblk );
    z_el     = reshape( GCOORD(3,EL2NOD(1:4,il:iu,:)), nnodel, nelblk );
    M_blk(:) = 0;

    for ip=1:nip % integration loop
        % Components of Jacobi matrices of all elements in block
        Jx       = x_el'*dNds{ip};
        Jy       = y_el'*dNds{ip};
        Jz       = z_el'*dNds{ip};
        % Determinants of Jacobi matrices ("Js")
        detJ     = Jx(:,1).*Jy(:,2).*Jz(:,3) ...
                 + Jx(:,2).*Jy(:,3).*Jz(:,1) ...
                 + Jx(:,3).*Jy(:,1).*Jz(:,2) ...
                 - Jx(:,3).*Jy(:,2).*Jz(:,1) ...
                 - Jx(:,1).*Jy(:,3).*Jz(:,2) ...
                 - Jx(:,2).*Jy(:,1).*Jz(:,3);
        if any(detJ<0)
            error('negative Jacobi')
        end

        weight = detJ*IP_w(ip); % Integration weight times area of triangle
        indx  = 1;
        for i=1:nnodel
            for j=i:nnodel
                M_blk(:,indx) = M_blk(:,indx) + weight .* (N{ip}(i)*N{ip}(j));
                indx = indx + 1;
            end
        end
    end % end of integration loop

    % Store element matrices for later assembly
    M_all(:,il:iu) = M_blk';
    
    il = il+nelblk;
    if (iblk==nblk-1)
        % Account for different number of elements in last block
        nelblk  = nel-iu; % number of remaining elements
        M_blk   = zeros(nelblk, nnodel*(nnodel+1)/2);
    end
    iu = iu + nelblk;
end % end of block loop

if exist('sparse_create','file')
    opts_mutils.symmetric  = 1;
    opts_mutils.n_node_dof = 1;
    opts_mutils.nthreads   = 1;
    MM = sparse_create(EL2NOD,M_all,opts_mutils);
else
    indx_j    = repmat(1:nnodel,nnodel,1); indx_i = indx_j';
    indx_i    = tril(indx_i); indx_i = indx_i(:); indx_i = indx_i(indx_i>0);
    indx_j    = tril(indx_j); indx_j = indx_j(:); indx_j = indx_j(indx_j>0);
    M_i       = EL2NOD(indx_i,:); M_i = M_i(:);
    M_j       = EL2NOD(indx_j,:); M_j = M_j(:);
    indx      = M_i < M_j;
    tmp       = M_j(indx);
    M_j(indx) = M_i(indx);
    M_i(indx) = tmp;
    MM        = sparse3(M_i, M_j, M_all);
end

MM = MM + tril(MM,-1)'; % other half is given by symmetry

end % END OF SUBFUNCTION assembly_mass_matrix_opt

% #########################################################################

function NV = assembly_rhs_opt(GCOORD,EL2NOD,N,dNds,IP_w,V_ip,nelblk)

nnodel = size(EL2NOD,1);
nel    = size(EL2NOD,2);
nelblk = min(nel,nelblk);
nblk   = ceil(nel/nelblk);
nip    = length(IP_w);

NV_blk = zeros(nelblk,nnodel);
NV_all = zeros(nnodel,nel);

il = 1;
iu = nelblk;
for iblk=1:nblk % loop over element blocks
    x_el      = reshape( GCOORD(1,EL2NOD(1:4,il:iu,:)), nnodel, nelblk );
    y_el      = reshape( GCOORD(2,EL2NOD(1:4,il:iu,:)), nnodel, nelblk );
    z_el      = reshape( GCOORD(3,EL2NOD(1:4,il:iu,:)), nnodel, nelblk );
    NV_blk(:) = 0;

    for ip=1:nip % integration loop
        % Components of Jacobi matrices of all elements in block
        Jx       = x_el'*dNds{ip};
        Jy       = y_el'*dNds{ip};
        Jz       = z_el'*dNds{ip};
        % Determinants of Jacobi matrices ("Js")
        detJ     = Jx(:,1).*Jy(:,2).*Jz(:,3) ...
                 + Jx(:,2).*Jy(:,3).*Jz(:,1) ...
                 + Jx(:,3).*Jy(:,1).*Jz(:,2) ...
                 - Jx(:,3).*Jy(:,2).*Jz(:,1) ...
                 - Jx(:,1).*Jy(:,3).*Jz(:,2) ...
                 - Jx(:,2).*Jy(:,1).*Jz(:,3);
        if any(detJ<0)
            error('negative Jacobi')
        end

        weight = detJ*IP_w(ip); % Integration weight times area of triangle
        NV_blk = NV_blk + (weight.*V_ip(il:iu,ip)) * N{ip}';
    end % end of integration loop

    % Store element matrices for later assembly
    NV_all(:,il:iu) = NV_blk';
    
    il = il+nelblk;
    if (iblk==nblk-1)
        % Account for different number of elements in last block
        nelblk   = nel-iu; % number of remaining elements
        NV_blk   = zeros(nelblk, nnodel);
    end
    iu = iu + nelblk;
end % end of block loop

NV  = accumarray(EL2NOD(:),NV_all(:));

end % END OF SUBFUNCTION assembly_rhs_opt