function ASSEMBLY = assembly_thermal_els_in_cone_iso_std...
    (MESH,EL2NOD,PhaseID,T,dt,VAR,num_scaling,ASSEMBLY,SF,PHYSICS,OPTS_T)

%===================================================================================================================================================================
% MODEL PARAMETERS
%===================================================================================================================================================================
GCOORD          = MESH.GCOORD;
els_in_cone_iso = MESH.els_in_cone_iso;
nnodel          = size(EL2NOD,1);
CC_all          = ASSEMBLY.CC_all;
i_all           = ASSEMBLY.i_all;
j_all           = ASSEMBLY.j_all;
rhs             = ASSEMBLY.rhs;
w_ip            = SF.w_ip;
N               = SF.N;
dNds            = SF.dNds;
RR_X_90_CCW     = [ 1  0  0 ; ...
                    0  0 -1 ; ...
                    0  1  0 ];

%===================================================================================================================================================================
% ELEMENT LOOP - MATRIX COMPUTATION
%===================================================================================================================================================================
for i = 1:size(els_in_cone_iso,2)
    iel = els_in_cone_iso(i);
    
    % Compute coordinates of each element inside cone and isoparametric in the 90° X-rotated spherical frame 
    GCOORD_this_el         = GCOORD(:,EL2NOD(:,iel));
    GCOORD_rot_this_el     = RR_X_90_CCW * GCOORD_this_el;
    GCOORD_SPH_rot_this_el = cartesian2spherical(GCOORD_rot_this_el);
    
    % Initialize element arrays
    CC_el  = zeros(nnodel);   % element coefficient matrix
    rhs_el = zeros(nnodel,1); % element r.h.s-vector
    elnods = EL2NOD(:,iel);

    %===============================================================================================================================================================
    % INTEGRATION LOOP
    %===============================================================================================================================================================
    for ip=1:OPTS_T.nip
        %===========================================================================================================================================================
        % CALCULATE 1st JACOBIAN (FROM SPHERICAL TO LOCAL COORDINATES) FOR ELEMENTS INSIDE THE CONE AND ISOPARAMETRIC
        % NOTE: These elements have curved edges in the 90° X-rotated spherical frame so the 1st Jacobian needs to be computed 
        %       at each integration point, i.e., calculated inside the integration loop. Quadratic 10-node shape functions 
        %       are needed to calculate the Jacobian since elements have curved edges in the 90° X-rotated spherical frame. 
        %===========================================================================================================================================================
        ECOORD  = GCOORD_SPH_rot_this_el; % 10-node element coordinates
        J_SL    = dNds{ip}'*ECOORD';      % Jacobian to transform spherical into local derivatives
        detJ_SL =  J_SL(1)*J_SL(5)*J_SL(9) + J_SL(4)*J_SL(8)*J_SL(3) + J_SL(7)*J_SL(2)*J_SL(6) ...
                 - J_SL(7)*J_SL(5)*J_SL(3) - J_SL(1)*J_SL(8)*J_SL(6) - J_SL(4)*J_SL(2)*J_SL(9);
        if detJ_SL<=0
            error('Error. \nNegative jacobian in element %d',iel);
        end
        
        %===========================================================================================================================================================
        % CALCULATE 2nd JACOBIAN (FROM CARTESIAN TO SPHERICAL COORDINATES) FOR ELEMENTS INSIDE THE CONE AND ISOPARAMETRIC 
        % NOTE: These elements have curved edges in the Cartesian frame so the 2nd Jacobian needs to be computed at each  
        %       integration point, i.e., calculated inside the integration loop. Quadratic 10-node shape functions are 
        %       needed to calculate the Jacobian since elements have curved edges in the 90° X-rotated spherical frame. 
        %===========================================================================================================================================================
        TH_gp    = sum(ECOORD.*repmat(N{ip},1,3)',2); % global coordinates (theta,phi,r) of the ip-th integration point 
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
%         J_CS(2,2)    = r_ip.*s_th.*c_ph;
%         J_CS(2,3)    =  0;
%         J_CS(3,1)    =  s_th.*c_ph;
%         J_CS(3,2)    =  s_th.*s_ph;
%         J_CS(3,3)    =  c_th;
%         detJ_CS_a    = det(J_CS);
%         diff_detJ_CS = detJ_CS - detJ_CS_a;
%         %---------------------------------------------------------------------------------------------------------------------------------------------------------
        
        %===========================================================================================================================================================
        % CALCULATE DOUBLE JACOBIAN
        %===========================================================================================================================================================
        % J_double is equal to J_SC * J_LS, where J_LS = inv(J_SL), so it would be J_SC * inv(J_SL), but it is faster and more accurate to use J_SC/J_SL
        % If we subtract J_SC*J_LS - J_SC/J_SL the result is the order of 1^-18 (machine precision) 
        J_double = J_SC/J_SL;
        
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
        
        %===========================================================================================================================================================
        % DERIVATIVES wrt GLOBAL COORDINATES
        %===========================================================================================================================================================
        dNdx     = J_double * dNds{ip}';      % global derivative of shape function
        
        %===========================================================================================================================================================
        % NUMERICAL INTEGRATION OF ELEMENT MATRICES
        %===========================================================================================================================================================
        weight   = detJ_CS*detJ_SL* w_ip(ip); % element area * Gauss weight (weights includes volume of master tetrahedron A=1/6)(V=0.5*A*h;A=area of base;h=height)
        NNt      = N{ip}*N{ip}'; % Ni x Nj
        if nnodel==4 && strcmp(OPTS_T.lumping,'yes')
            % Lumped mass matrix
            NNt = diag(sum(NNt,2)); % Lumped
        end
        
        %===========================================================================================================================================================
        % PROPERTIES OF ELEMENTS AT ip-TH EVALUATION POINT
        %===========================================================================================================================================================
        Cond_el  = calc_element_conductivity ...
            (VAR,OPTS_T,PHYSICS,EL2NOD,PhaseID,iel,N{ip});
        RhoCp_el = calc_element_RhoCp...
            (VAR,OPTS_T,PHYSICS,EL2NOD,PhaseID,iel,N{ip});
        dQdt_el  = calc_element_dQdt...
            (VAR,OPTS_T,PHYSICS,EL2NOD,PhaseID,iel,N{ip});
        
        %===========================================================================================================================================================
        % CALCULATE ELEMENT MATRICES
        %===========================================================================================================================================================
        CC_el  = CC_el  + (dNdx'*dNdx)*num_scaling*dt*Cond_el*weight + NNt*weight*RhoCp_el;
        rhs_el = rhs_el + NNt*weight*T(elnods)*RhoCp_el;
        
        % SOURCE TERM
        if ~isempty(dQdt_el)
            if num_scaling~=1
                error('Non-zero source term must be verified if num_scaling~=1 !!');
            end
            rhs_el = rhs_el + N{ip}*(num_scaling*dt*dQdt_el.*weight)'; % FIXME now scaling with num_scaling to make similar to thermal2d
        end
    end % END OF INTEGRATION LOOP
    
    %======================================================================
    % RECOVER MATRIX FORM APPLYING ROTATIONS
    %======================================================================
    % Recover the matrix CC and vector rhs from rotated counterparts  
    % because the elements have been rotated 90° around X axis
    % NOTE: rhs vecotres are the same since they are independent of the 
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
    
    %===============================================================================================================================================================
    % STORE DATA OF ALL ELEMENTS FOR ASSEMBLY
    %===============================================================================================================================================================
    rows          = double(elnods)*ones(1,nnodel);  % row location in global matrices
    cols          = ones(nnodel,1)*double(elnods)'; % column location in global matrices
    CC_all(:,iel) = CC_el(:);
    i_all(:,iel)  = rows(:); 
    j_all(:,iel)  = cols(:);
    
    %======================================================================
    % R.H.S.-VECTOR CAN BE ACCUMULATED DIRECTLY
    %======================================================================
    rhs(elnods)  = rhs(elnods) + rhs_el;
end % END OF ELEMENT LOOP

ASSEMBLY.CC_all = CC_all;
ASSEMBLY.i_all  = i_all;
ASSEMBLY.j_all  = j_all;
ASSEMBLY.rhs    = rhs;

end % END OF FUNCTION assembly_thermal_els_in_cone_iso_std