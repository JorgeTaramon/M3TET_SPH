function ASSEMBLY = assembly_stokes_els_out_cone_cross_2pi_opt...
    (ASSEMBLY,SF,VAR,MESH,PHYSICS,NUMSCALE,OPTS)

% MILAMIN method (blocks of elements handled at once)
% see Dabrowski et al. 2008 (doi:10.1029/2007GC001719)

%==========================================================================
% MODEL PARAMETERS
%==========================================================================
els_selct  = MESH.els_out_cone_cross_2pi;
nel        = length(els_selct);
EL2NOD     = uint32(MESH.EL2NOD{1}(:  ,els_selct)); % connectivity matrix for velocity problem
EL2NODP    = uint32(MESH.EL2NOD{1}(1:4,els_selct)); % connectivity matrix for pressure problem
GCOORD_SPH = MESH.GCOORD_SPH; % node coordinates
els_ref    = MESH.els_ref;    % elements in refined region
ndim       = 3;               % number of spatial dimensions
PhaseID    = MESH.PhaseID(els_selct);
if length(PhaseID)==1; PhaseID = ones(nel,1); end
nnodel     = size(EL2NOD,1);  % number of nodes in each element
nvertx     = MESH.nvertx;     % number of vertices (or corners) in an element
nUdofel    = ndim*nnodel;     % number of velocity dofs in each element
nPdofel    = size(EL2NODP,1); % number of pressure dofs in each element
IP_w       = SF.IP_w;
NU         = SF.NU;
dNUds      = SF.dNUds;
NP         = SF.NP;
dN4ds      = SF.dN4ds;

%==========================================================================
% STORAGE FOR DATA OF ALL ELEMENT MATRICES/VECTORS
%==========================================================================
Fb_all  = zeros(nUdofel              , nel);  % storage for global force vector
if isfield(VAR,'DilEl')
    Fp_all  = zeros(nPdofel          , nel);  % storage for global dilation force vector
end
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
il        = 1;
iu        = nelblk;

tic
%==========================================================================
% BLOCK LOOP - MATRIX COMPUTATION
%==========================================================================
for ib=1:nblk % loop over element blocks
    els_blk = els_selct(il:iu);
    
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
        Dil_blk = VAR.DilEl(els_blk); % dilation prescribed for each element
    else
        Dil_blk = [];
    end
    
    %================================================================================================================================
    % CALCULATE 1st JACOBIAN (FROM SPHERICAL TO LOCAL COORDINATES) FOR ELEMENTS OUTSIDE THE CONE AND CROSSING PHI = 2PI 
    % NOTE: These elements have straight edges in the 180° Z-rotated spherical frame so the 1st Jacobian is the same 
    %       for all integration points, i.e., calculated once before the integration loop. Further, linear 4-node shape
    %       functions are sufficient to calculate the Jacobian. 
    %================================================================================================================================
    VCOORD_th = reshape( GCOORD_SPH(1,EL2NOD(1:nvertx,il:iu)), nvertx, nelblk );
    VCOORD_ph = reshape( GCOORD_SPH(2,EL2NOD(1:nvertx,il:iu)), nvertx, nelblk );
    VCOORD_r  = reshape( GCOORD_SPH(3,EL2NOD(1:nvertx,il:iu)), nvertx, nelblk );
    
    % 180° counterclocwise rotation around z-axis (adding pi to phi)
    VCOORD_ph      = VCOORD_ph + pi;
    ind            = VCOORD_ph>2*pi;
    VCOORD_ph(ind) = VCOORD_ph(ind) - 2*pi;
    
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
    for ip=1:OPTS.nip
        %============================================================================================================================
        % CALCULATE 2nd JACOBIAN (FROM CARTESIAN TO SPHERICAL COORDINATES) FOR ELEMENTS OUTSIDE THE CONE AND CROSSING PHI = 2PI 
        % NOTE: These elements have curved edges in the Cartesian frame so the 2nd Jacobian needs to be computed at each  
        %       integration point, i.e., calculated inside the integration loop. Linear 4-node shape functions are sufficient 
        %       to calculate the Jacobian since elements have straight edges in the 180° Z-rotated spherical frame. 
        %============================================================================================================================
        th_ip    = VCOORD_th'*NP{ip};
        ph_ip    = VCOORD_ph'*NP{ip};
        r_ip     = VCOORD_r' *NP{ip};
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
        % CALCULATE DOUBLE JACOBIAN AND STIFFNESS MATRICES
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
        dNUdx  = Jx_dbl*dNUds{ip}';
        dNUdy  = Jy_dbl*dNUds{ip}';
        dNUdz  = Jz_dbl*dNUds{ip}';
        
        weight = detJ_CS.*detJ_SL.*IP_w(ip); % integration weight times volumes of tetrahedrons

        %==================================================================
        % PROPERTIES OF ELEMENTS AT ip-TH INTEGRATION POINT
        %==================================================================
        Dens_blk = calc_element_dens(VAR,OPTS,PHYSICS,EL2NOD,PhaseID,il:iu,NP{ip},th_ip,ph_ip,r_ip);
        Visc_blk = calc_element_visc(VAR,OPTS,PHYSICS,EL2NOD,PhaseID,il:iu,NP{ip},th_ip,ph_ip,r_ip);
        
        if OPTS.return_dens_el
            DensEl(il:iu,ip) = Dens_blk(:);
        end
        if OPTS.return_visc_ip
            ViscIP(il:iu,ip) = Visc_blk(:);
        end
        
        % Gravitational force at ip-th integration point
        % Note that Bscale is new defined: Bscale = L0^2/(Visc0 * U0) 
        % (it does no longer includes grav. acceleration and density)
%         Fg_blk = PHYSICS.g * NUMSCALE.Bscale .* (Dens_blk-NUMSCALE.Dens0);
%         if ~isempty(els_ref) % mesh with embedded high resolution subregion
%             Fg_blk(~ismember(els_blk,els_ref)) = 0; % compute gravitational force only in the refined region
%                                                     % (Fg outside of the ref region is equal to zero)
%         end
        
        if ~isempty(els_ref)
            I = ismember(els_blk,els_ref); % compute gravitational force only in the refined region
            J = ismember(els_ref,els_blk);
            Fg_blk    = zeros(size(Dens_blk,1),1); % gravitational force out of the refined region is zero
            Fg_blk(I) = PHYSICS.g * NUMSCALE.Bscale * (Dens_blk(I) - NUMSCALE.Dens0_T_bary(J)); 
        else
            Fg_blk = PHYSICS.g * NUMSCALE.Bscale .* (Dens_blk-NUMSCALE.Dens0);
        end
        
        %==================================================================
        % NUMERICAL INTEGRATION OF ELEMENT MATRICES
        %==================================================================
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
            tmp1        = weight.*NP_blk(:,i);
            tmp2        = tmp1(:,ones(1,nnodel));
            ii          = (i-1)*nUdofel + (1:3:nUdofel-2);
            G_blk(:,ii) = G_blk(:,ii) + tmp2.*dNUdx;
            ii          = (i-1)*nUdofel + (2:3:nUdofel-1);
            G_blk(:,ii) = G_blk(:,ii) + tmp2.*dNUdy;
            ii          = (i-1)*nUdofel + (3:3:nUdofel);
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
        Fb_blk(:,1:3:nUdofel) = Fb_blk(:,1:3:nUdofel) + (s_th.*c_ph .* weight .* Fg_blk) * NU{ip}';
        Fb_blk(:,2:3:nUdofel) = Fb_blk(:,2:3:nUdofel) + (s_th.*s_ph .* weight .* Fg_blk) * NU{ip}';
        Fb_blk(:,3:3:nUdofel) = Fb_blk(:,3:3:nUdofel) + (c_th       .* weight .* Fg_blk) * NU{ip}';
        
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
    K_all(:,il:iu)  = K_blk';
    G_all(:,il:iu)  = G_blk';
    M_all(:,il:iu)  = M_blk';
    Fb_all(:,il:iu) = Fb_blk';
    if ~isempty(Dil_blk)
        Fp_all(:,il:iu)  = Fp_blk';
    end
    
    %======================================================================
    % READJUST START, END AND SIZE OF BLOCK
    %======================================================================
    il = il + nelblk;
    if ib==nblk-1
        nelblk   = nel-iu;
    end
    iu  = iu  + nelblk; % update counter for elements in block
end % END OF BLOCK LOOP

%======================================================================
% RECOVER MATRIX FORM APPLYING ROTATIONS
%======================================================================
% Recover the matrices KK, GG and Fb from rotated counterparts because
% the elements have been rotated 180° around Z axis
% NOTE: MM matrices are the same since they are independent of the
%       element coordinates

% KK matrix
% For each element: KK_unrot = [RR_Z_180_CCW]' * KK * [RR_Z_180_CCW]
% 
% The pattern for the 3 degrees of freedom (dofs) for each node is 
% (compare against standard assembly):
%                  1  2  3 | 4  ... (column dofs)
%
%            1     1  1 -1 | 1
%   row dofs 2     1  1 -1 | 1
%            3    -1 -1  1 |-1
%            --------------|------
%            4     1  1 -1 | 1
%

% Generate a full matrix that for marks where each element's stiffness
% matrix has to change sign:
Ae                  = ones(nUdofel);
Ae(3:3:nUdofel,:)   = -Ae(3:3:nUdofel,:);
Ae(:,3:3:nUdofel)   = -Ae(:,3:3:nUdofel);
% Now extract the lower triangular part because only this part has been
% calculated above
indx_tril           = reshape(1:nUdofel^2,nUdofel,nUdofel);
indx_tril           = tril(indx_tril);
indx_tril           = indx_tril(indx_tril>0);
Ae                  = Ae(indx_tril);
indx_minus          = find(Ae==-1); % <--- these are the rows in K_all that
                                    %      have to change sign
K_all(indx_minus,:) = -K_all(indx_minus,:);

% GG matrix
G_all(1:3:nUdofel*nPdofel,:) = -G_all(1:3:nUdofel*nPdofel,:);
G_all(2:3:nUdofel*nPdofel,:) = -G_all(2:3:nUdofel*nPdofel,:);

% Fb vector
Fb_all(1:3:nUdofel,:) = -Fb_all(1:3:nUdofel,:);
Fb_all(2:3:nUdofel,:) = -Fb_all(2:3:nUdofel,:);

% STORE VALUES FOR LATER ASSEMBLY
ASSEMBLY.KKv  (:,els_selct) = K_all;
ASSEMBLY.GGv  (:,els_selct) = G_all;
ASSEMBLY.MMv  (:,els_selct) = M_all;
ASSEMBLY.Fb   (:,els_selct) = Fb_all;
if isfield(VAR,'DilEl')
    ASSEMBLY.Fpdil(:,els_selct) = Fp_all;
end

% Average density in each element is sufficient
if OPTS.return_dens_el
    ASSEMBLY.DensEl(els_selct,:) = DensEl;
end
if OPTS.return_visc_ip
    ASSEMBLY.ViscEl(els_selct,:) = ViscEl;
end

end % END OF FUNCTION assembly_stokes_els_out_cone_cross_2pi_opt