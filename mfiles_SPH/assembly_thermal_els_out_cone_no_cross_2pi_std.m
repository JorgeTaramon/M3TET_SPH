function ASSEMBLY = assembly_thermal_els_out_cone_no_cross_2pi_std...
    (MESH,EL2NOD,PhaseID,T,dt,VAR,num_scaling,ASSEMBLY,SF,PHYSICS,OPTS_T)

%===================================================================================================================================================================
% MODEL PARAMETERS
%===================================================================================================================================================================
GCOORD_SPH                = MESH.GCOORD_SPH;
els_out_cone_no_cross_2pi = MESH.els_out_cone_no_cross_2pi;
nnodel                    = size(EL2NOD,1);
CC_all                    = ASSEMBLY.CC_all;
i_all                     = ASSEMBLY.i_all;
j_all                     = ASSEMBLY.j_all;
rhs                       = ASSEMBLY.rhs;
w_ip                      = SF.w_ip;
N                         = SF.N;
dNds                      = SF.dNds;
N4                        = SF.N4;
dN4ds                     = SF.dN4ds;

%===================================================================================================================================================================
% ELEMENT LOOP - MATRIX COMPUTATION
%===================================================================================================================================================================
for i = 1:size(els_out_cone_no_cross_2pi,2)
    iel = els_out_cone_no_cross_2pi(i);
    %===============================================================================================================================================================
    % CALCULATE 1st JACOBIAN (FROM SPHERICAL TO LOCAL COORDINATES) FOR ELEMENTS OUTSIDE THE CONE AND NOT CROSSING PHI = 2PI 
    % NOTE: These elements have straight edges in the original spherical frame so the 1st Jacobian is the same 
    %       for all integration points, i.e., calculated once before the integration loop. Further, linear 4-node 
    %       shape functions are sufficient to calculate the Jacobian. 
    %===============================================================================================================================================================
    VCOORD  = GCOORD_SPH(:,EL2NOD(1:4,iel)); % element node vertex coordinates
    J_SL    = dN4ds*VCOORD'; % Jacobian to transform spherical into local derivatives
    detJ_SL =  J_SL(1)*J_SL(5)*J_SL(9) + J_SL(4)*J_SL(8)*J_SL(3) + J_SL(7)*J_SL(2)*J_SL(6) ...
             - J_SL(7)*J_SL(5)*J_SL(3) - J_SL(1)*J_SL(8)*J_SL(6) - J_SL(4)*J_SL(2)*J_SL(9);
    if detJ_SL<=0
        error('Error. \nNegative jacobian in element %d',iel);
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
    CC_el  = zeros(nnodel);   % element coefficient matrix
    rhs_el = zeros(nnodel,1); % element r.h.s-vector
    elnods = EL2NOD(:,iel);

    %===============================================================================================================================================================
    % INTEGRATION LOOP
    %===============================================================================================================================================================
    for ip=1:OPTS_T.nip
        %===========================================================================================================================================================
        % CALCULATE 2nd JACOBIAN (FROM CARTESIAN TO SPHERICAL COORDINATES) FOR ELEMENTS OUTSIDE THE CONE AND NOT CROSSING PHI = 2PI 
        % NOTE: These elements have curved edges in the Cartesian frame so the 2nd Jacobian needs to be computed at each 
        %       integration point, i.e., calculated inside the integration loop. Linear 4-node shape functions are sufficient 
        %       to calculate the Jacobian since elements have straight edges in the original spherical frame. 
        %===========================================================================================================================================================
        TH_gp    = sum(VCOORD.*repmat(N4{ip},1,3)',2); % global coordinates (theta,phi,r) of the ip-th integration point 
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

end % END OF FUNCTION assembly_thermal_els_out_cone_no_cross_2pi_std