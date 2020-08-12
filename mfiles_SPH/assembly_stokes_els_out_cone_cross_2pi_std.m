function ASSEMBLY = assembly_stokes_els_out_cone_cross_2pi_std...
    (ASSEMBLY,SF,D,B,EL2DOF,VAR,MESH,PHYSICS,NUMSCALE,OPTS)

%==========================================================================
% MODEL PARAMETERS
%==========================================================================
GCOORD                 = MESH.GCOORD;
PhaseID                = MESH.PhaseID;
EL2NOD                 = double(MESH.EL2NOD{1});
EL2NODP                = double(MESH.EL2NOD{1}(1:4,:));
els_ref                = MESH.els_ref;
ndim                   = 3;
nnodel                 = size(EL2NOD,1);
nUdofel                = ndim*nnodel;
nPdofel                = size(EL2NODP,1);
els_out_cone_cross_2pi = MESH.els_out_cone_cross_2pi;
KKv                    = ASSEMBLY.KKv;
KKi                    = ASSEMBLY.KKi;
KKj                    = ASSEMBLY.KKj;
Fb                     = ASSEMBLY.Fb;
GGv                    = ASSEMBLY.GGv;
GGi                    = ASSEMBLY.GGi;
GGj                    = ASSEMBLY.GGj;
MMv                    = ASSEMBLY.MMv;
MMi                    = ASSEMBLY.MMi;
MMj                    = ASSEMBLY.MMj;
Fpdil                  = ASSEMBLY.Fpdil;
DensEl                 = ASSEMBLY.DensEl;
ViscIP                 = ASSEMBLY.ViscIP;
IP_w                   = SF.IP_w;
NU                     = SF.NU;
dNUds                  = SF.dNUds;
NP                     = SF.NP;
dN4ds                  = SF.dN4ds;
RR_Z_180_CCW           = [-1  0  0 ; ...
                           0 -1  0 ; ...
                           0  0  1 ];

%==========================================================================
% ELEMENT LOOP - MATRIX COMPUTATION
%==========================================================================
for i = 1:size(els_out_cone_cross_2pi,2)
    iel = els_out_cone_cross_2pi(i);
    %===============================================================================================================================================================
    % CALCULATE 1st JACOBIAN (FROM SPHERICAL TO LOCAL COORDINATES) FOR ELEMENTS OUTSIDE THE CONE AND CROSSING PHI = 2PI 
    % NOTE: These elements have straight edges in the 180° Z-rotated spherical frame so the 1st Jacobian is the same 
    %       for all integration points, i.e., calculated once before the integration loop. Further, linear 4-node shape
    %       functions are sufficient to calculate the Jacobian. 
    %===============================================================================================================================================================
    GCOORD_this_el         = GCOORD(:,EL2NOD(:,iel));
    GCOORD_rot_this_el     = RR_Z_180_CCW * GCOORD_this_el;
    GCOORD_SPH_rot_this_el = cartesian2spherical(GCOORD_rot_this_el);
    VCOORD                 = GCOORD_SPH_rot_this_el(:,1:4); % element node vertex coordinates
    J_SL                   = dN4ds*VCOORD'; % Jacobian to transform spherical into local derivatives
    detJ_SL                =  J_SL(1)*J_SL(5)*J_SL(9) + J_SL(4)*J_SL(8)*J_SL(3) + J_SL(7)*J_SL(2)*J_SL(6) ...
                            - J_SL(7)*J_SL(5)*J_SL(3) - J_SL(1)*J_SL(8)*J_SL(6) - J_SL(4)*J_SL(2)*J_SL(9);
    if detJ_SL<=0
        error('Error. \nNegative jacobian in element %s',iel);
    end
    
%     % Alternative way to compute J_LS (useful for debugging) -----------------------------------------------------------------------------------------------------
%     
%     Th14      = VCOORD(1,1) - VCOORD(1,4);
%     Th24      = VCOORD(1,2) - VCOORD(1,4);
%     Th34      = VCOORD(1,3) - VCOORD(1,4);
%     Ph14      = VCOORD(2,1) - VCOORD(2,4);
%     Ph24      = VCOORD(2,2) - VCOORD(2,4);
%     Ph34      = VCOORD(2,3) - VCOORD(2,4);
%     R14       = VCOORD(3,1) - VCOORD(3,4);
%     R24       = VCOORD(3,2) - VCOORD(3,4);
%     R34       = VCOORD(3,3) - VCOORD(3,4);
%     detJ_SL_a = Th14*(Ph24*R34 - Ph34*R24) + ...
%                 Th24*(Ph34*R14 - Ph14*R34) + ...
%                 Th34*(Ph14*R24 - Ph24*R14);
%     J_LS_a    = [Ph24*R34  - Ph34*R24    Ph34*R14  - Ph14*R34    Ph14*R24  - Ph24*R14; ...
%                  Th34*R24  - Th24*R34    Th14*R34  - Th34*R14    Th24*R14  - Th14*R24; ...
%                  Th24*Ph34 - Th34*Ph24   Th34*Ph14 - Th14*Ph34   Th14*Ph24 - Th24*Ph14] / detJ_SL_a;
%     % J_LS_a is the Jacobian to transform local into spherical derivatives
%     % J_LS_a is equal to (J_SL)^-1
%     I_1st_double_Jac = J_SL*J_LS_a; % it should be the identity
%     weight_double    = zeros(1,OPTS.nip);
%     %-------------------------------------------------------------------------------------------------------------------------------------------------------------
    
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
        %==================================================================
        % PROPERTIES OF ELEMENTS AT ip-TH INTEGRATION POINT
        %==================================================================
        Dens_iel = calc_element_dens(VAR,OPTS,PHYSICS,EL2NOD,PhaseID,iel,NP{ip});
        Visc_iel = calc_element_visc(VAR,OPTS,PHYSICS,EL2NOD,PhaseID,iel,NP{ip});
        
        if OPTS.return_dens_el
            DensEl(iel,ip) = Dens_iel(:);
        end
        if OPTS.return_visc_ip
            ViscIP(iel,ip) = Visc_iel(:);
        end
        
        % Gravitational force at ip-th integration point. Note that Bscale is new defined: Bscale = L0^2/(Visc0 * U0) 
        % It does no longer include grav. acceleration and density
        if ~isempty(els_ref) % mesh with embedded high resolution subregion
            [I,J] = ismember(iel,els_ref);
            if I % compute gravitational force only in the refined region
%                 Fg_el = PHYSICS.g * NUMSCALE.Bscale * (Dens_iel-NUMSCALE.Dens0);
                
                Fg_el = PHYSICS.g * NUMSCALE.Bscale * (Dens_iel - NUMSCALE.Dens0_T_bary(J));
            else
                Fg_el = 0; % gravitational force out of the refined region is zero
            end
        else % regular mesh -> compute gravitational force at each element
            Fg_el = PHYSICS.g * NUMSCALE.Bscale * (Dens_iel-NUMSCALE.Dens0);
        end
        
        %===========================================================================================================================================================
        % CALCULATE 2nd JACOBIAN (FROM CARTESIAN TO SPHERICAL COORDINATES) FOR ELEMENTS OUTSIDE THE CONE AND CROSSING PHI = 2PI 
        % NOTE: These elements have curved edges in the Cartesian frame so the 2nd Jacobian needs to be computed at each  
        %       integration point, i.e., calculated inside the integration loop. Linear 4-node shape functions are sufficient 
        %       to calculate the Jacobian since elements have straight edges in the 180° Z-rotated spherical frame. 
        %===========================================================================================================================================================
        TH_gp    = sum(VCOORD.*repmat(NP{ip},1,3)',2); % global coordinates (theta,phi,r) of the ip-th integration point 
        th_ip    = TH_gp(1); % theta evaluated at ip-th integration point
        ph_ip    = TH_gp(2); % phi evaluated at ip-th integration point
        r_ip     = TH_gp(3); % r evaluated at ip-th integration point
        if th_ip == 0
            error('Error. \nip %s in element %s has theta = 0 (it makes the det(J_CS) = Inf)',num2str(ip),num2str(iel))
        end
        s_th     = sin(th_ip);
        c_th     = cos(th_ip);
        s_ph     = sin(ph_ip);
        c_ph     = cos(ph_ip);
        
        detJ_CS  = r_ip.^2.*s_th;
        J_SC     = [r_ip.*s_th.*c_th.*c_ph   -r_ip.*s_ph   r_ip.^2.*s_th.^2.*c_ph; ...
                    r_ip.*s_th.*c_th.*s_ph    r_ip.*c_ph   r_ip.^2.*s_th.^2.*s_ph; ...
                            -r_ip.*s_th.^2        0           r_ip.^2.*s_th.*c_th] / detJ_CS;
        
%         % Alternative way to compute J_CS (useful for debugging) -------------------------------------------------------------------------------------------------
%         
%         % J_SC is the jacobian to transform spherical into Cartesian derivatives
%         % J_SC is equal to (J_CS)^-1
%         J_CS         = zeros(3);
%         J_CS(1,1)    =  r_ip.*c_th.*c_ph;
%         J_CS(1,2)    =  r_ip.*c_th.*s_ph;
%         J_CS(1,3)    = -r_ip.*s_th;
%         J_CS(2,1)    = -r_ip.*s_th.*s_ph;
%         J_CS(2,2)    =  r_ip.*s_th.*c_ph;
%         J_CS(2,3)    =  0;
%         J_CS(3,1)    =  s_th.*c_ph;
%         J_CS(3,2)    =  s_th.*s_ph;
%         J_CS(3,3)    =  c_th;
%         detJ_CS_a    = det(J_CS);
%         diff_detJ_CS = detJ_CS - detJ_CS_a;
%         %---------------------------------------------------------------------------------------------------------------------------------------------------------
        
        %===========================================================================================================================================================
        % CALCULATE DOUBLE JACOBIAN AND STIFFNESS MATRICES
        %===========================================================================================================================================================
        % J_double is equal to J_SC * J_LS, where J_LS = inv(J_SL), so it would be J_SC * inv(J_SL), but it is faster and more accurate to use J_SC/J_SL
        % If we subtract J_SC*J_LS - J_SC/J_SL the result is the order of 1^-18 (machine precision) 
        J_double = J_SC/J_SL;
        dNUdx    = J_double * dNUds{ip}';     % global derivative of shape function
        weight   = detJ_CS*detJ_SL* IP_w(ip); % element area * Gauss weight (weights includes volume of master tetrahedron A=1/6)(V=0.5*A*h;A=area of base;h=height)
        
%         % Alternative way to compute J_double (useful for debugging) ---------------------------------------------------------------------------------------------
%         
%         J_double_a      = zeros(3);
%         J_double_a(1,1) =  ( ((c_th.*c_ph)./r_ip).*(Ph24*R34  - Ph34*R24) + (-s_ph./(r_ip.*s_th)).*(Th34*R24  - Th24*R34) + s_th.*c_ph.*(Th24*Ph34 - Th34*Ph24) );
%         J_double_a(1,2) =  ( ((c_th.*c_ph)./r_ip).*(Ph34*R14  - Ph14*R34) + (-s_ph./(r_ip.*s_th)).*(Th14*R34  - Th34*R14) + s_th.*c_ph.*(Th34*Ph14 - Th14*Ph34) );
%         J_double_a(1,3) =  ( ((c_th.*c_ph)./r_ip).*(Ph14*R24  - Ph24*R14) + (-s_ph./(r_ip.*s_th)).*(Th24*R14  - Th14*R24) + s_th.*c_ph.*(Th14*Ph24 - Th24*Ph14) );
%         J_double_a(2,1) =  ( ((c_th.*s_ph)./r_ip).*(Ph24*R34  - Ph34*R24) + ( c_ph./(r_ip.*s_th)).*(Th34*R24  - Th24*R34) + s_th.*s_ph.*(Th24*Ph34 - Th34*Ph24) );
%         J_double_a(2,2) =  ( ((c_th.*s_ph)./r_ip).*(Ph34*R14  - Ph14*R34) + ( c_ph./(r_ip.*s_th)).*(Th14*R34  - Th34*R14) + s_th.*s_ph.*(Th34*Ph14 - Th14*Ph34) );
%         J_double_a(2,3) =  ( ((c_th.*s_ph)./r_ip).*(Ph14*R24  - Ph24*R14) + ( c_ph./(r_ip.*s_th)).*(Th24*R14  - Th14*R24) + s_th.*s_ph.*(Th14*Ph24 - Th24*Ph14) );
%         J_double_a(3,1) =  (        (-s_th./r_ip).*(Ph24*R34  - Ph34*R24) +                       0                       +       c_th.*(Th24*Ph34 - Th34*Ph24) );
%         J_double_a(3,2) =  (        (-s_th./r_ip).*(Ph34*R14  - Ph14*R34) +                       0                       +       c_th.*(Th34*Ph14 - Th14*Ph34) );
%         J_double_a(3,3) =  (        (-s_th./r_ip).*(Ph14*R24  - Ph24*R14) +                       0                       +       c_th.*(Th14*Ph24 - Th24*Ph14) );
%         
%         J_double_a        = J_double_a/detJ_SL_a;
%         diff_double_Jac   = J_double-J_double_a;
%         weight_double(ip) = weight;
%         %---------------------------------------------------------------------------------------------------------------------------------------------------------
        
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
        rhatGP           = [s_th*c_ph; s_th*s_ph; c_th]; % unit vector in radial direction (we don not use [-s_th*c_ph; -s_th*s_ph; -c_th] 
                                                         % because g is already negative (-9.81))
        Fbe(1:3:nUdofel) = Fbe(1:3:nUdofel) + NU{ip} .* Fg_el * rhatGP(1,:)' * weight; % buoyancy force x-dir
        Fbe(2:3:nUdofel) = Fbe(2:3:nUdofel) + NU{ip} .* Fg_el * rhatGP(2,:)' * weight; % buoyancy force y-dir
        Fbe(3:3:nUdofel) = Fbe(3:3:nUdofel) + NU{ip} .* Fg_el * rhatGP(3,:)' * weight; % buoyancy force z-dir
        
%         % Alternative way to compute force vector (useful for debugging) -----------------------------------------------------------------------------------------
%         
%         ECOORD = GCOORD(:,EL2NOD(:,iel));
%         xyz_ip  = ECOORD*NU{ip};          % physical coordinates of integration point
%         r_ip    = sqrt(sum(xyz_ip.^2,1)); % radius of integration point
%         rhat_ip = xyz_ip ./ r_ip;         % scale them all by l
%         rhatGP - rhat_ip
%             
%         VCOORD_debug     = GCOORD(:,EL2NOD(1:4,iel)); % element node vertex coordinates
%         X_gp             = sum(VCOORD_debug.*repmat(NP{ip},1,3)',2); % global coordinates (x,y,z) of the ip-th integration point
%         scale            = sum(X_gp.^2,1);   % scale = radius^2
%         rhatGP_debug     = X_gp ./ repmat( sqrt(scale) ,[3 1]); % scale them all by l
%         Fbe_debug        = zeros(size(Fbe));
%         Fbe_debug(1:3:nUdofel) = Fbe_debug(1:3:nUdofel) + NU{ip} .* Fg_el * rhatGP_debug(1,:)' * weight; % buoyancy force x-dir
%         Fbe_debug(2:3:nUdofel) = Fbe_debug(2:3:nUdofel) + NU{ip} .* Fg_el * rhatGP_debug(2,:)' * weight; % buoyancy force y-dir
%         Fbe_debug(3:3:nUdofel) = Fbe_debug(3:3:nUdofel) + NU{ip} .* Fg_el * rhatGP_debug(3,:)' * weight; % buoyancy force z-dir
%         Fbe - Fbe_debug
%         %---------------------------------------------------------------------------------------------------------------------------------------------------------
        
        % R.h.s. term on incompressibility equation from non-zero dilation
        Fpe              = Fpe + NP{ip}*DilEl*weight;
    end % END OF INTEGRATION LOOP
    
    %======================================================================
    % RECOVER MATRIX FORM APPLYING ROTATIONS
    %======================================================================
    % Recover the matrices KK, GG and Fb from KK_rot, GG_rot and Fb_rot for those elements that have been rotated 180° around Z axis
    % NOTE: MM and MM_rot are the same since they are independent of the element coordinates

    % KK matrix
    % For each element: KK_unrot = [RR_Z_180_CCW]' * KK * [RR_Z_180_CCW]
    KKe_temp = KKe;
    ndim     = 3;
    nnodel   = size(EL2NOD,1);
    for j = 1:3:ndim*nnodel
        KKe_temp(1:3:end,j)   =  KKe(1:3:end,j);   % (1,1) is equal to  (1,1)
        KKe_temp(2:3:end,j)   =  KKe(2:3:end,j);   % (2,1) is equal to  (2,1)
        KKe_temp(3:3:end,j)   = -KKe(3:3:end,j);   % (3,1) is equal to -(3,1)
        KKe_temp(1:3:end,j+1) =  KKe(1:3:end,j+1); % (1,2) is equal to  (1,2)
        KKe_temp(2:3:end,j+1) =  KKe(2:3:end,j+1); % (2,2) is equal to  (2,2)
        KKe_temp(3:3:end,j+1) = -KKe(3:3:end,j+1); % (3,2) is equal to -(3,2)
        KKe_temp(1:3:end,j+2) = -KKe(1:3:end,j+2); % (1,3) is equal to -(1,3)
        KKe_temp(2:3:end,j+2) = -KKe(2:3:end,j+2); % (2,3) is equal to -(2,3)
        KKe_temp(3:3:end,j+2) =  KKe(3:3:end,j+2); % (3,3) is equal to  (3,3)
    end
    KKe = KKe_temp;

    % GG matrix
    GGe_temp            =  GGe;
    GGe_temp(1:3:end,:) = -GGe(1:3:end,:);
    GGe_temp(2:3:end,:) = -GGe(2:3:end,:);
    GGe                 =  GGe_temp;

    % Fb vector
    Fbe_temp            =  Fbe;
    Fbe_temp(1:3:end,:) = -Fbe(1:3:end,:);
    Fbe_temp(2:3:end,:) = -Fbe(2:3:end,:);
    Fbe                 =  Fbe_temp;
    
    %======================================================================
    % STORE DATA OF ALL ELEMENTS IN BLOCK FOR ASSEMBLY
    %======================================================================
    gUdof = EL2DOF(:,iel);  % global velocity dofs of current element
    gPdof = EL2NODP(:,iel); % global pressure dofs of current element
    
    % Store stiffness matrix of each element for assembly
    rows         = gUdof*ones(1,nUdofel);
    cols         = ones(nUdofel,1)*gUdof';
    KKv(:,iel)   = KKe(:);  % values
    KKi(:,iel)   = rows(:); % i indices
    KKj(:,iel)   = cols(:); % j indices

    % Assemble global force vector including body force term
    Fb(gUdof)    = Fb(gUdof) + Fbe;
    
    % Assemble global pressure-dilation r.h.s. vector
    Fpdil(gPdof) = Fpdil(gPdof) + Fpe;
    
    % Store gradient matrix of each element for assembly
    rows         = gUdof*ones(1,nPdofel);  % rows of GGe are velocity dofs
    cols         = ones(nUdofel,1)*gPdof'; % cols of GGe are pressure dofs 
    GGv(:,iel)   = GGe(:);  % values
    GGi(:,iel)   = rows(:); % i indices
    GGj(:,iel)   = cols(:); % j indices
    
    % Store (1/viscosity)-scaled mass matrix of each element
    rows         = gPdof*ones(1,nPdofel);
    cols         = ones(nPdofel,1)*gPdof';
    MMv(:,iel)   = MMe(:);  % values
    MMi(:,iel)   = rows(:); % i indices
    MMj(:,iel)   = cols(:); % j indices
    
end % END OF ELEMENT LOOP

% Average density in each element is sufficient
if OPTS.return_dens_el
    DensEl = mean(DensEl,2);
end

ASSEMBLY.KKv    = KKv;
ASSEMBLY.KKi    = KKi;
ASSEMBLY.KKj    = KKj;
ASSEMBLY.Fb     = Fb;
ASSEMBLY.GGv    = GGv;
ASSEMBLY.GGi    = GGi;
ASSEMBLY.GGj    = GGj;
ASSEMBLY.MMv    = MMv;
ASSEMBLY.MMi    = MMi;
ASSEMBLY.MMj    = MMj;
ASSEMBLY.Fpdil  = Fpdil;
ASSEMBLY.DensEl = DensEl;
ASSEMBLY.ViscIP = ViscIP;

end % END OF FUNCTION assembly_stokes_els_out_cone_cross_2pi_std