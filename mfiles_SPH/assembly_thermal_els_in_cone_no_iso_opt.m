function ASSEMBLY = assembly_thermal_els_in_cone_no_iso_opt...
    (MESH,EL2NOD,PhaseID,T,dt,VAR,num_scaling,ASSEMBLY,SF,PHYSICS,OPTS_T)

% MILAMIN method (blocks of elements handled at once)
% see Dabrowski et al. 2008 (doi:10.1029/2007GC001719)

%==========================================================================
% MODEL PARAMETERS
%==========================================================================
els_selct  = MESH.els_in_cone_no_iso;
nel        = length(els_selct);
EL2NOD     = uint32(EL2NOD(:,els_selct));
GCOORD     = MESH.GCOORD;
nnodel     = size(EL2NOD,1);
nvertx     = MESH.nvertx;
w_ip       = SF.w_ip;
N          = SF.N;
dNds       = SF.dNds;
N4         = SF.N4;
dN4ds      = SF.dN4ds;

%==========================================================================
% STORAGE FOR DATA OF ALL ELEMENT MATRICES/VECTORS
%==========================================================================
CC_all     = zeros(nnodel*(nnodel+1)/2,nel);
rhs_all    = zeros(nnodel,nel);

%==========================================================================
% BLOCKING PARAMETERS (nelblk must be < nel)
%==========================================================================
nelblk    = min(nel, OPTS_T.nelblk); % in case nel<nelblk
nblk      = ceil(nel/nelblk);        % number of blocks
il        = 1;
iu        = nelblk;

%==========================================================================
% BLOCK LOOP - MATRIX COMPUTATION
%==========================================================================
for ib=1:nblk % loop over element blocks
    
    %======================================================================
    % STORAGE FOR DATA OF ELEMENTS IN BLOCK
    %======================================================================
    CC_blk  = zeros(nelblk, nnodel*(nnodel+1)/2);
    rhs_blk = zeros(nelblk, nnodel);
    T_blk   = T(EL2NOD(:,il:iu));
    
    %===============================================================================================================================================================
    % CALCULATE 1st JACOBIAN (FROM SPHERICAL TO LOCAL COORDINATES) FOR ELEMENTS INSIDE THE CONE AND NOT ISOPARAMETRIC 
    % NOTE: These elements have straight edges in the 90� X-rotated spherical frame so the 1st Jacobian is the same 
    %       for all integration points, i.e., calculated once before the integration loop. Further, linear 4-node shape
    %       functions are sufficient to calculate the Jacobian. 
    %===============================================================================================================================================================
    % Rotation matrix:
    % RR_X_90_CCW  = [ 1  0  0 ;
    %                  0  0 -1 ;
    %                  0  1  0 ];
    % ==> ECOORD_x = -ECOORD_x;
    % ==> ECOORD_y = -ECOORD_z;
    % ==> ECOORD_z =  ECOORD_y;
    VCOORD_x_rot = reshape( GCOORD(1,EL2NOD(1:nvertx,il:iu)), nvertx, nelblk );
    VCOORD_y_rot = reshape(-GCOORD(3,EL2NOD(1:nvertx,il:iu)), nvertx, nelblk );
    VCOORD_z_rot = reshape( GCOORD(2,EL2NOD(1:nvertx,il:iu)), nvertx, nelblk );

    [VCOORD_th,VCOORD_ph,VCOORD_r] = elblk_cartesian2spherical...
        (VCOORD_x_rot,VCOORD_y_rot,VCOORD_z_rot);
    
    Jth_SL    = VCOORD_th'*dN4ds;
    Jph_SL    = VCOORD_ph'*dN4ds;
    Jr_SL     = VCOORD_r' *dN4ds;
    % J_SL of 1st element in block would be:
    %       [Jth_SL(1,:)' Jph_SL(1,:)' Jr_SL(1,:)']
    
    detJ_SL   = elblk_detA(Jth_SL,Jph_SL,Jr_SL);
    if any(detJ_SL<0)
        error('negative Jacobi')
    end
    
    % Inversion of Jacobi matrices "J_SL" of all elements in block
    % (J_LS = inv(J_SL))
    [Jth_LS,Jph_LS,Jr_LS] = elblk_invert_A(Jth_SL,Jph_SL,Jr_SL,detJ_SL);
    % J_LS (=inv(J_SL) of 1st element in block would be:
    %       [Jth_LS(1,:); Jph_LS(1,:); Jr_LS(1,:)]
    
    %======================================================================
    % INTEGRATION LOOP: DEAL WITH A BLOCK OF ELEMENTS AT ONCE
    %======================================================================
    for ip=1:OPTS_T.nip
        %===========================================================================================================================================================
        % CALCULATE 2nd JACOBIAN (FROM CARTESIAN TO SPHERICAL COORDINATES) FOR ELEMENTS INSIDE THE CONE AND NOT ISOPARAMETRIC 
        % NOTE: These elements have curved edges in the Cartesian frame so the 2nd Jacobian needs to be computed at each integration 
        %       point, i.e., calculated inside the integration loop. Linear 4-node shape functions are sufficient to calculate 
        %       the Jacobian since elements have straight edges in the 90� X-rotated spherical frame. 
        %===========================================================================================================================================================
        th_ip    = VCOORD_th'*N4{ip};
        ph_ip    = VCOORD_ph'*N4{ip};
        r_ip     = VCOORD_r' *N4{ip};
        if any(th_ip == 0)
            error('Error. \nip %i in element %i has theta = 0 (it makes the det(J_CS) = Inf)',ip,els_blk(th_ip==0))
        end
        s_th     = sin(th_ip);
        c_th     = cos(th_ip);
        s_ph     = sin(ph_ip);
        c_ph     = cos(ph_ip);
        
        detJ_CS  = r_ip.^2.*s_th;
        
        J1_SC    = [ r_ip.*s_th.*c_th.*c_ph  -r_ip.*s_ph      r_ip.^2.*s_th.^2.*c_ph] ./ repmat(detJ_CS,1,3);
        J2_SC    = [ r_ip.*s_th.*c_th.*s_ph   r_ip.*c_ph      r_ip.^2.*s_th.^2.*s_ph] ./ repmat(detJ_CS,1,3);
        J3_SC    = [-r_ip.*s_th.^2           zeros(nelblk,1)  r_ip.^2.*s_th   .*c_th] ./ repmat(detJ_CS,1,3);
        % J_SC of 1st element in block would be: [J1_SC(1,:);
        %                                         J2_SC(1,:);
        %                                         J3_SC(1,:)]
        
        %============================================================================================================================
        % CALCULATE DOUBLE JACOBIAN
        %============================================================================================================================
        % J_double is equal to J_SC * J_LS, where J_LS = inv(J_SL)
        % i.e.: J_double = J_SC * inv(J_SL);
        [Jx_dbl,Jy_dbl,Jz_dbl] = elblk_AtxBt(J1_SC,J2_SC,J3_SC,Jth_LS,Jph_LS,Jr_LS);
        % J_double of 1st element in block would be:
        %       [Jx_dbl(1,:)' Jy_dbl(1,:)' Jz_dbl(1,:)']
        
        [Jx_dbl,Jy_dbl,Jz_dbl] = elblk_transpose_A(Jx_dbl,Jy_dbl,Jz_dbl);
        % Now J_double of 1st element in block would be: [Jx_dbl(1,:);
        %                                                 Jy_dbl(1,:);
        %                                                 Jz_dbl(1,:)]
        
        %==================================================================
        % DERIVATIVES wrt GLOBAL COORDINATES
        %==================================================================
        dNdx  = Jx_dbl*dNds{ip}';
        dNdy  = Jy_dbl*dNds{ip}';
        dNdz  = Jz_dbl*dNds{ip}';
        
        %==========================================================
        % NUMERICAL INTEGRATION OF ELEMENT MATRICES
        %==========================================================
        weight = detJ_CS.*detJ_SL.*w_ip(ip); % integration weight times volumes of tetrahedrons
        Ni     = N{ip};
        NNt    = Ni*Ni'; % Ni x Nj
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
                CC_blk(:,indx) = CC_blk(:,indx) ...
                   + (dNdx(:,i).*dNdx(:,j) + dNdy(:,i).*dNdy(:,j) + dNdz(:,i).*dNdz(:,j)) ...
                      .* num_scaling .* dt .* Cond_blk ...
                   + NNt(i,j).*RhoCp_blk;
                indx = indx + 1;
            end
        end
% %         % To show full CC-matrix of first 4-node element in block:
% %         CC_el = CC_blk(1,:);CC_el([1 2 3 4; 2 5 6 7; 3 6 8 9; 4 7 9 10])
        
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
    end % END OF INTEGRATION LOOP
    
    %==============================================================
    % WRITE DATA INTO GLOBAL STORAGE
    %==============================================================
    CC_all(:,il:iu)	 = CC_blk';
    rhs_all(:,il:iu) = rhs_blk';
    
    %==============================================================
    % READJUST START, END AND SIZE OF BLOCK. REALLOCATE MEMORY
    %==============================================================
    il  = il + nelblk;
    if(ib==nblk-1)
        nelblk = nel-iu;
    end
    iu  = iu + nelblk;
end % END OF BLOCK LOOP

%======================================================================
% RECOVER MATRIX FORM APPLYING ROTATIONS
%======================================================================
% Recover the matrix CC and vector rhs from rotated counterparts
% because the elements have been rotated 90� around X axis
% NOTE: rhs vectors are the same since they are independent of the
%       element coordinates
%       CC matrices are the same due to the way they are computed:
%               CC  = (dNdx'*dNdx) * kappa * weight
%       where "kappa" and "weight" are independent of the element
%       coodinates. The global derivaties are given by:
%                   dNdx = J_double * dNds{ip}'
%       where dNds{ip} is independent of element coordinates and
%       J_double is dependent of the element coordinates. This means
%       that dNdx is dependent of the element coordinates (i.e., dNdx
%       will be different in each frame: original and rotated).
%       However, when computing dNdx'*dNdx, the result is independent
%       of the frame.

% STORE VALUES FOR LATER ASSEMBLY
ASSEMBLY.CC_v  (:,els_selct) = CC_all;
ASSEMBLY.rhs_v (:,els_selct) = rhs_all;
ASSEMBLY.nblk = ASSEMBLY.nblk + nblk;

end % END OF FUNCTION assembly_thermal_els_in_cone_no_iso_opt