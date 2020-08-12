function [CC,MM] = assembly_diffusion_3d(GCOORD,EL2NOD,kappa,OPTS_D,fidl)
% Usage: [CC,MM] = assembly_diffusion_3d(GCOORD,EL2NOD,kappa,OPTS_D,fidl)
%
% Purpose: 
%   Calculate the element matrices and vectors; assemble global
%   counterparts 
%
% Input:
%   GCOORD   : [matrix]    : Cartesian coordinates of all nodes in mesh
%   EL2NOD   : [matrix]    : connectivity matrix
%   kappa    : [scalar]    : diffusivity (must be in correct units !!!)
%   OPTS_D   : [structure] : options for assembly procedure
%   fidl     : [scalar]    : 
%
% Output:
%   CC       : [sparsemat] : global conductivity (diffusivity) matrix
%   MM       : [sparsemat] : global mass (capacity) matrix
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)

% JH March 2011
% JH March 2013
% JH & JMT Nov 2016 : cleaned up, restructured

switch OPTS_D.method
    case 'std' % standard element assembly: good for understanding FE
        [CC,MM] = assembly_std(GCOORD,EL2NOD,kappa,OPTS_D,fidl); % *SUBFUNCTION*
    case 'opt' % a faster block-wise assembly (MILAMIN style)
        [CC,MM] = assembly_opt(GCOORD,EL2NOD,kappa,OPTS_D,fidl); % *SUBFUNCTION*
% %         [CC2,MM2] = assembly_std(GCOORD,EL2NOD,kappa,OPTS_D,fidl); % *SUBFUNCTION*
% %         disp(full(max(max(abs(CC-CC2))) / max(max((abs(CC))))))
% %         disp(full(max(max(abs(MM-MM2))) / max(max((abs(MM))))))
end

end % END OF FUNCTION assembly_diffusion_3d

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function [CC,MM] = assembly_std(GCOORD,EL2NOD,kappa,OPTS_D,fidl)

t      = tic;

%==========================================================================
% MODEL INFO
%==========================================================================
nel    = size(EL2NOD,2);
nnodel = size(EL2NOD,1);

%==========================================================================
% STORAGE FOR GLOBAL MATRICES AND VECTORS
%==========================================================================
C_all  = zeros(nnodel*nnodel,nel); % storage for conductivity matrix coefficients
M_all  = zeros(nnodel*nnodel,nel); % storage for capacity matrix coefficients
i_all  = zeros(nnodel*nnodel,nel); % storage for the matrix i indixes
j_all  = zeros(nnodel*nnodel,nel); % storage for the matrix j indixes

%==========================================================================
% PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
%==========================================================================
[x_ip,w_ip] = ip_tetrahedron(OPTS_D.nip);
[N,dNds]    = sf_dsf_tet(x_ip,nnodel,'cell');

%==========================================================================
% ELEMENT LOOP - MATRIX COMPUTATION
%==========================================================================
for iel = 1:nel
    % Initialize element arrays
    CC_el = zeros(nnodel); % element conductivity matrix
    MM_el = zeros(nnodel); % element capacity matrix
    
    %==============================================================
    % INTEGRATION LOOP
    %==============================================================
    for ip=1:OPTS_D.nip
        %==========================================================
        % LOAD SHAPE FUNCTIONS DERIVATIVES FOR INTEGRATION POINT
        %==========================================================
        Ni    = N{ip};
        dNids = dNds{ip};
            
        %==================================================================
        % CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
        %==================================================================
        elnods = EL2NOD(:,iel);
        ECOORD = GCOORD(:,elnods); % element node coordinates
        J      = dNids'*ECOORD';   % Jacobi matrix (mapping reference element ==> current element)
        detJ   =   J(1)*J(5)*J(9) + J(4)*J(8)*J(3) + J(7)*J(2)*J(6) ...
                 - J(7)*J(5)*J(3) - J(1)*J(8)*J(6) - J(4)*J(2)*J(9); % Determinate of Jacobi matrix
        if(detJ<=0)
            error('Negative jacobian in element %1i',iel);
        end
        
        %==========================================================
        % DERIVATIVES wrt GLOBAL COORDINATES
        %==========================================================
        dNdx   = J\dNds{ip}';   % global derivatives of SF 
        
        %==========================================================
        % NUMERICAL INTEGRATION OF ELEMENT MATRICES
        %==========================================================
        weight = w_ip(ip)*detJ;
        
        % calculate element matrices
        CC_el = CC_el + (dNdx'*dNdx) * kappa * weight;
        MM_el = MM_el + (Ni*Ni') * weight;
        
    end % END OF INTEGRATION LOOP
    
    %======================================================================
    % STORE DATA OF ALL ELEMENTS FOR ASSEMBLY
    %======================================================================
    rows         = double(elnods)*ones(1,nnodel);  % row location in global matrices
    cols         = ones(nnodel,1)*double(elnods)'; % column location in global matrices
    C_all(:,iel) = CC_el(:);
    M_all(:,iel) = MM_el(:);
    i_all(:,iel) = rows(:); 
    j_all(:,iel) = cols(:);
end % END OF ELEMENT LOOP
t_elloop = toc(t);t=tic;

%==========================================================================
% MERGE ELEMENT MATRICES INTO GLOBAL SPARSE MATRIX
%==========================================================================
CC  = sparse3(i_all(:),j_all(:),C_all(:)); % global conductivity matrix
MM  = sparse3(i_all(:),j_all(:),M_all(:)); % global conductivity matrix
t_sparse = toc(t);

if fidl
    fprintf(fidl,'   element loop                          : %7.2f sec \n',t_elloop);
    fprintf(fidl,'   making sparse matrix                  : %7.2f sec\n',t_sparse);
end

end % END OF SUBFUNCTION assembly_std

% #########################################################################

function [CC,MM] = assembly_opt(GCOORD,EL2NOD,kappa,OPTS_D,fidl)

t       = tic;

%==========================================================================
% MODEL INFO
%==========================================================================
nel     = size(EL2NOD,2);
nnodel  = size(EL2NOD,1);

%==========================================================================
% PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
%==========================================================================
[x_ip,w_ip] = ip_tetrahedron(OPTS_D.nip);
[N,dNds]    = sf_dsf_tet(x_ip,nnodel,'cell');

%==========================================================================
% DECLARE VARIABLES (ALLOCATE MEMORY)
%==========================================================================
C_all   = zeros(nnodel*(nnodel+1)/2,nel);
M_all   = zeros(nnodel*(nnodel+1)/2,nel);

%==========================================================================
% BLOCKING PARAMETERS (nelblk must be < nel)
%==========================================================================
nelblk   = min(nel, OPTS_D.nelblk);
nblk     = ceil(nel/nelblk);
il       = 1;
iu       = nelblk;

%==================================================================
% BLOCK LOOP - MATRIX COMPUTATION
%==================================================================
for ib = 1:nblk
    C_blk = zeros(nelblk, nnodel*(nnodel+1)/2);
    M_blk = zeros(nelblk, nnodel*(nnodel+1)/2);
    
    %==============================================================
    % INTEGRATION LOOP
    %==============================================================
    for ip=1:OPTS_D.nip
        %==========================================================
        % LOAD SHAPE FUNCTIONS DERIVATIVES FOR INTEGRATION POINT
        %==========================================================
        Ni          = N{ip};
        dNids       = dNds{ip};
            
        %======================================================================
        % CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
        %======================================================================
        [detJ,invJx,invJy,invJz] = calc_jacobian(GCOORD,EL2NOD(:,il:iu),dNds{ip});
        
        %==========================================================
        % DERIVATIVES wrt GLOBAL COORDINATES
        %==========================================================
        dNdx      = invJx*dNids';
        dNdy      = invJy*dNids';
        dNdz      = invJz*dNids';
        
        %==========================================================
        % NUMERICAL INTEGRATION OF ELEMENT MATRICES
        %==========================================================
        weight    = w_ip(ip)*detJ;
        
        %==============================================================
        % PROPERTIES OF ELEMENTS AT ip-TH EVALUATION POINT
        %==============================================================
        indx = 1;
        for i = 1:nnodel
            for j = i:nnodel
                % Conductivity matrices
                C_blk(:,indx) = C_blk(:,indx) + kappa .* weight .* ...
                                                (dNdx(:,i).*dNdx(:,j) + ...
                                                 dNdy(:,i).*dNdy(:,j) + ...
                                                 dNdz(:,i).*dNdz(:,j));
                % Capacity matrices
                M_blk(:,indx) = M_blk(:,indx) + weight .* Ni(i) * Ni(j);
                indx = indx + 1;
            end
        end
% %         % To show full C-matrix of first 4-node element in block
% %         C_el = C_blk(1,:);C_el([1 2 3 4; 2 5 6 7; 3 6 8 9; 4 7 9 10])
    end
    
    %==============================================================
    % WRITE DATA INTO GLOBAL STORAGE
    %==============================================================
    C_all(:,il:iu) = C_blk';
    M_all(:,il:iu) = M_blk';
    
    %==============================================================
    % READJUST START, END AND SIZE OF BLOCK. REALLOCATE MEMORY
    %==============================================================
    il  = il + nelblk;
    if(ib==nblk-1)
        nelblk = nel-iu;
    end
    iu  = iu + nelblk;
end
t_elloop = toc(t);t=tic;

%==========================================================================
% MERGE ELEMENT MATRICES INTO GLOBAL SPARSE MATRIX
%==========================================================================
opts_mutils.symmetric  = 1;
opts_mutils.n_node_dof = 1;
opts_mutils.nthreads   = 1; % >1 is very slow (????)
CC = sparse_create(EL2NOD,C_all,opts_mutils);
MM = sparse_create(EL2NOD,M_all,opts_mutils);
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