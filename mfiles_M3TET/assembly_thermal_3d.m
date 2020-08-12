function [CC,rhs] = assembly_thermal_3d(GCOORD,EL2NOD,PhaseID,T,dt,VAR,num_scaling,PHYSICS,OPTS_T,fidl)
% Usage: [CC,rhs] = assembly_thermal_3d(GCOORD,EL2NOD,PhaseID,T,dt,VAR,num_scaling,PHYSICS,OPTS_T,fidl)
% 
% Purpose: 
%   Calculate the element matrices and vectors; assemble global
%   counterparts 
%
% Input:
%   GCOORD      : [matrix]    : Cartesian coordinates of all nodes in mesh
%   EL2NOD      : [matrix]    : connectivity matrix
%   PhaseID     : [vector]    : phase indices
%   T           : [colvector] : variable field to be diffused
%   dt          : [scalar]    : time over which is diffused
%   VAR         : [structure] : major variable fields
%   num_scaling : [scalar]    : scaling factor to macth the non-SI units
%   PHYSICS     : [structure] : physical properties
%   OPTS_D      : [structure] : options for assembly procedure
%   fidl        : [scalar]    : 
% 
% Output:
%   CC          : [matrix]    : global conductivity matrix
%   rhs         : [vector]    : right hand side vector
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Aug 2013
% JH & JMT Nov 2016 : cleaned up, restructured

switch OPTS_T.method
    case 'std'
        [CC,rhs] = assembly_std(GCOORD,EL2NOD,PhaseID,T,dt,VAR,num_scaling,PHYSICS,OPTS_T,fidl); % *SUBFUNCTION*
    case 'opt'
        [CC,rhs] = assembly_opt(GCOORD,EL2NOD,PhaseID,T,dt,VAR,num_scaling,PHYSICS,OPTS_T,fidl); % *SUBFUNCTION*
% %         [CC2,rhs2] = assembly_std(GCOORD,EL2NOD,PhaseID,T,dt,VAR,num_scaling,PHYSICS,OPTS_T,fidl); % *SUBFUNCTION*
% %         disp(full(max(max(abs(CC-CC2))) / max(max((abs(CC))))))
% %         disp(max(abs(rhs-rhs2)) / max(abs(rhs)))
end

end % END OF FUNCTION assembly_thermal_3d

% #########################################################################
%                               SUBFUNCTIONS
% #########################################################################

function [CC,rhs] = assembly_std(GCOORD,EL2NOD,PhaseID,T,dt,VAR,num_scaling,PHYSICS,OPTS_T,fidl)

t      = tic;

%==========================================================================
% MODEL INFO
%==========================================================================
nel    = size(EL2NOD,2);
nnodel = size(EL2NOD,1);
nnod   = length(T);

%==========================================================================
% STORAGE FOR GLOBAL MATRICES AND VECTORS
%==========================================================================
C_all  = zeros(nnodel*nnodel,nel); % storage for conductivity matrix coefficients
i_all  = zeros(nnodel*nnodel,nel); % storage for the matrix i indixes
j_all  = zeros(nnodel*nnodel,nel); % storage for the matrix j indixes
rhs    = zeros(nnod,1);            % storage for rhs force vector coefficients

%==========================================================================
% PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
%==========================================================================
[x_ip,w_ip] = ip_tetrahedron(OPTS_T.nip);
[N,dNds]    = sf_dsf_tet(x_ip,nnodel,'cell');

%==========================================================================
% ELEMENT LOOP - MATRIX COMPUTATION
%==========================================================================
for iel = 1:nel
    % Initialize element arrays
    CC_el  = zeros(nnodel);   % element coefficient matrix
    rhs_el = zeros(nnodel,1); % element r.h.s-vector
    
    %==============================================================
    % INTEGRATION LOOP
    %==============================================================
    for ip=1:OPTS_T.nip
        %==========================================================
        % LOAD SHAPE FUNCTIONS DERIVATIVES FOR INTEGRATION POINT
        %==========================================================
        Ni          = N{ip};
        dNids       = dNds{ip};
            
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
        NNt    = Ni*Ni'; % Ni x Nj
        if nnodel==4 && strcmp(OPTS_T.lumping,'yes')
            % Lumped mass matrix
            NNt = diag(sum(NNt,2)); % Lumped
        end
        
        %==============================================================
        % PROPERTIES OF ELEMENTS AT ip-TH EVALUATION POINT
        %==============================================================
        Cond_el  = calc_element_conductivity ...
            (VAR,OPTS_T,PHYSICS,EL2NOD,PhaseID,iel,Ni);
        RhoCp_el = calc_element_RhoCp...
            (VAR,OPTS_T,PHYSICS,EL2NOD,PhaseID,iel,Ni);
        dQdt_el  = calc_element_dQdt...
            (VAR,OPTS_T,PHYSICS,EL2NOD,PhaseID,iel,Ni);

        % calculate element matrices
        CC_el  = CC_el  + (dNdx'*dNdx)*num_scaling*dt*Cond_el*weight + NNt*weight*RhoCp_el;
        rhs_el = rhs_el + NNt*weight*T(elnods)*RhoCp_el;
        
        % SOURCE TERM
        if ~isempty(dQdt_el)
            if num_scaling~=1
                error('Non-zero source term must be verified if num_scaling~=1 !!');
            end
            rhs_el = rhs_el + Ni*(num_scaling*dt*dQdt_el.*weight)'; % FIXME now scaling with num_scaling to make similar to thermal2d
        end
    end % END OF INTEGRATION LOOP
    
    %======================================================================
    % STORE DATA OF ALL ELEMENTS FOR ASSEMBLY
    %======================================================================
    rows         = double(elnods)*ones(1,nnodel);  % row location in global matrices
    cols         = ones(nnodel,1)*double(elnods)'; % column location in global matrices
    C_all(:,iel) = CC_el(:);
    i_all(:,iel) = rows(:); 
    j_all(:,iel) = cols(:);
    
    %======================================================================
    % R.H.S.-VECTOR CAN BE ACCUMULATED DIRECTLY
    %======================================================================
    rhs(elnods)  = rhs(elnods) + rhs_el;
end % END OF ELEMENT LOOP
t_elloop = toc(t);t=tic;

%==========================================================================
% MERGE ELEMENT MATRICES INTO GLOBAL SPARSE MATRIX
%==========================================================================
CC  = sparse3(i_all(:),j_all(:),C_all(:)); % global conductivity matrix
CC  = tril(CC); % this is to have the same CC than in 'opt' method
t_sparse = toc(t);

if fidl
    fprintf(fidl,'   element loop                          : %7.2f sec \n',t_elloop);
    fprintf(fidl,'   making sparse matrix                  : %7.2f sec\n',t_sparse);
end

end % END OF SUBFUNCTION assembly_std

% #########################################################################

function [CC,rhs] = assembly_opt(GCOORD,EL2NOD,PhaseID,T,dt,VAR,num_scaling,PHYSICS,OPTS_T,fidl)

t      = tic;

%==========================================================================
% MODEL INFO
%==========================================================================
nel     = size(EL2NOD,2);
nnodel  = size(EL2NOD,1);

%==========================================================================
% PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
%==========================================================================
[x_ip,w_ip] = ip_tetrahedron(OPTS_T.nip);
[N,dNds]    = sf_dsf_tet(x_ip,nnodel,'cell');

%==========================================================================
% DECLARE VARIABLES (ALLOCATE MEMORY)
%==========================================================================
C_all   = zeros(nnodel*(nnodel+1)/2,nel);
rhs_all = zeros(nnodel,nel);

%==========================================================================
% BLOCKING PARAMETERS (nelblk must be < nel)
%==========================================================================
nelblk   = min(nel, OPTS_T.nelblk);
nblk     = ceil(nel/nelblk);
il       = 1;
iu       = nelblk;

%==================================================================
% BLOCK LOOP - MATRIX COMPUTATION
%==================================================================
for ib = 1:nblk
    C_blk	 = zeros(nelblk, nnodel*(nnodel+1)/2);
    rhs_blk  = zeros(nelblk, nnodel);
    T_blk    = T(EL2NOD(:,il:iu));
    
    %==============================================================
    % INTEGRATION LOOP
    %==============================================================
    for ip=1:OPTS_T.nip
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
        NNt       = Ni*Ni'; % Ni x Nj
        if nnodel==4 && strcmp(OPTS_T.lumping,'yes')
            % Lumped mass matrix
            NNt = diag(sum(NNt,2)); % Lumped
        end
        
        %==============================================================
        % PROPERTIES OF ELEMENTS AT ip-TH EVALUATION POINT
        %==============================================================
        Cond_blk  = calc_element_conductivity ...
            (VAR,OPTS_T,PHYSICS,EL2NOD,PhaseID,il:iu,Ni);
        RhoCp_blk = calc_element_RhoCp...
            (VAR,OPTS_T,PHYSICS,EL2NOD,PhaseID,il:iu,Ni);
        dQdt_blk  = calc_element_dQdt...
            (VAR,OPTS_T,PHYSICS,EL2NOD,PhaseID,il:iu,Ni);

        % Cond_blk and RhoCp_blk now include the integration weights!
        Cond_blk  = weight .* Cond_blk;
        RhoCp_blk = weight .* RhoCp_blk;
        
        indx = 1;
        for i = 1:nnodel
            for j = i:nnodel
                % Note: Cond_blk and RhoCp_blk include integration weights!
                C_blk(:,indx) = C_blk(:,indx) ...
                   + (dNdx(:,i).*dNdx(:,j) + dNdy(:,i).*dNdy(:,j) + dNdz(:,i).*dNdz(:,j)) ...
                      .* num_scaling .* dt .* Cond_blk ...
                   + NNt(i,j).*RhoCp_blk;
                indx = indx + 1;
            end
        end
% %         % To show full C-matrix of first 4-node element in block:
% %         C_el = C_blk(1,:);C_el([1 2 3 4; 2 5 6 7; 3 6 8 9; 4 7 9 10])
        
        % RIGHT HAND SIDE
        % Note: Cond_blk and RhoCp_blk include integration weights!
        rhs_blk = rhs_blk + ( NNt*(T_blk.*(RhoCp_blk*ones(1,nnodel))' ))';
        
        % SOURCE TERM
        if ~isempty(dQdt_blk)
            if num_scaling~=1
                error('Non-zero source term must be verified if num_scaling~=1 !!');
            end
            rhs_blk = rhs_blk + ( Ni*(num_scaling*dt*dQdt_blk.*weight)' )'; % FIXME now scaling with num_scaling to make similar to thermal2d
        end
    end
    
    %==============================================================
    % WRITE DATA INTO GLOBAL STORAGE
    %==============================================================
    C_all(:,il:iu)	 = C_blk';
    rhs_all(:,il:iu) = rhs_blk';
    
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
% ASSEMBLE R.H.S. VECTOR
%==========================================================================
rhs = accumarray(EL2NOD(:), rhs_all(:));

%==========================================================================
% MERGE ELEMENT MATRICES INTO GLOBAL SPARSE MATRIX
%==========================================================================
opts_mutils.symmetric  = 1;
opts_mutils.n_node_dof = 1;
opts_mutils.nthreads   = 1; % >1 is very slow (????)
CC = sparse_create(EL2NOD,C_all,opts_mutils);
t_sparse = toc(t);

if fidl
    fprintf(fidl,'   element loop (nelblk=%4i,nblk=%6i): %7.2f sec \n',...
        OPTS_T.nelblk,nblk,t_elloop);
    fprintf(fidl,'   making sparse matrix                  : %7.2f sec\n',t_sparse);
end

end % END OF SUBFUNCTION assembly_opt