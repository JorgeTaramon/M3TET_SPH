function VAR = calc_stresses_3d(VAR,MESH,SETTINGS,PHYSICS,NUMSCALE,dt,nelblk)

% =========================================================================
% MODEL PARAMETERS
% =========================================================================
ndim    = 3;            % number of spatial dimensions
nel     = MESH.nel;     % number of elements
GCOORD  = MESH.GCOORD;
EL2NOD  = uint32(MESH.EL2NOD{1});  % connectivity matrix for velocity problem
PhaseID = MESH.PhaseID;
if length(PhaseID)==1; PhaseID = ones(MESH.nel,1); end
nnodel  = size(EL2NOD,1);  % number of nodes in each element
nvertx  = MESH.nvertx;     % number of vertices (or corners) in an element
SETTINGS.is_elastic = 'no';

%==========================================================================
% VELOCITY VECTOR
%==========================================================================
U          = zeros(3*MESH.nnod,1);
U(1:3:end) = VAR.Ux; % velocity vector Ux(1), Uy(1), Uz(1), Ux(2), Uy(2), ...
U(2:3:end) = VAR.Uy; % velocity vector Ux(1), Uy(1), Uz(1), Ux(2), Uy(2), ...
U(3:3:end) = VAR.Uz; % velocity vector Ux(1), Uy(1), Uz(1), Ux(2), Uy(2), ...

%==========================================================================
% PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
%==========================================================================
nip_stress = 4; % size(VAR.Stress_xx,2);
[IP_X,~]   = ip_tetrahedron(nip_stress);
  % local coordinates and weights of points for integration of
  % velocity/pressure matrices
[NU,dNUds] = sf_dsf_tet(IP_X,nnodel,'cell');
  % velocity shape functions and their derivatives
[NL,dNLds] = sf_dsf_tet(IP_X,nvertx,'cell');
  % derivatives of linear (4-node) shape functions are used to calculate
  % each element's Jacobian

  
%==========================================================================
% BLOCKING PARAMETERS (nelblk must be < nel)
%==========================================================================
nelblk    = min(nel,nelblk); % in case nel<nelblk
nblk      = ceil(nel/nelblk);     % number of blocks
il        = 1;
iu        = nelblk;


% =========================================================================
% BEGIN OF POST PROCESSING (STRESSES etc)
% =========================================================================
for iblk=1:nblk % LOOP OVER ELEMENT BLOCKS
    %======================================================================
    % CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
    %======================================================================
    % NOTE: For tetrahedral elements with non-curved edges the Jacobian is
    %       the same for all integration points so that it can be 
    %       calculated once before the integration loop. Linear 4-node 
    %       shape functions are sufficient to calculate the Jacobian.
    [~,invJx,invJy,invJz] = calc_jacobian(GCOORD,EL2NOD(1:nvertx,il:iu),dNLds{1});

    %======================================================================
    % STRESS CALCULATION LOOP
    %======================================================================
    eldof0  = ndim*EL2NOD(:,il:iu)-3;
    Ux_blk  = reshape(U(eldof0+1)',nelblk,nnodel);
    Uy_blk  = reshape(U(eldof0+2)',nelblk,nnodel);
    Uz_blk  = reshape(U(eldof0+3)',nelblk,nnodel);
    
    for ip=1:nip_stress % LOOP OVER STRESS EVALUATION POINTS
        %==================================================================
        % COORDINATES OF CORNER/VERTEX NODES OF ALL ELEMENTS IN BLOCK
        %==================================================================
        ECOORD_x = reshape( GCOORD(1,EL2NOD(:,il:iu)), nnodel, nelblk );
        ECOORD_y = reshape( GCOORD(2,EL2NOD(:,il:iu)), nnodel, nelblk );
        ECOORD_z = reshape( GCOORD(3,EL2NOD(:,il:iu)), nnodel, nelblk );
        
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
        
        %==================================================================
        % VISCOSITY OF ELEMENTS AT ip-TH EVALUATION POINT
        %==================================================================
        Visc_blk = calc_element_visc...
            (VAR,SETTINGS,PHYSICS,EL2NOD,PhaseID,il:iu,NL{ip});
        
        %==================================================================
        % DERIVATIVES wrt GLOBAL COORDINATES
        %==================================================================
        dNUdx   = invJx*dNUds{ip}';
        dNUdy   = invJy*dNUds{ip}';
        dNUdz   = invJz*dNUds{ip}';

        dUxdx   = sum(Ux_blk.*dNUdx,2);
        dUxdy   = sum(Ux_blk.*dNUdy,2);
        dUxdz   = sum(Ux_blk.*dNUdz,2);
        dUydx   = sum(Uy_blk.*dNUdx,2);
        dUydy   = sum(Uy_blk.*dNUdy,2);
        dUydz   = sum(Uy_blk.*dNUdz,2);
        dUzdx   = sum(Uz_blk.*dNUdx,2);
        dUzdy   = sum(Uz_blk.*dNUdy,2);
        dUzdz   = sum(Uz_blk.*dNUdz,2);

        Exy_blk = 0.5*(dUxdy + dUydx);
        Exz_blk = 0.5*(dUxdz + dUzdx);
        Eyz_blk = 0.5*(dUydz + dUzdy);
 
        %==================================================================
        % VISCOUS STRESSES
        %==================================================================
        Txx_blk =  4/3.*Visc_blk.*dUxdx - 2/3.*Visc_blk.*dUydy - 2/3.*Visc_blk.*dUzdz;
        Tyy_blk = -2/3.*Visc_blk.*dUxdx + 4/3.*Visc_blk.*dUydy - 2/3.*Visc_blk.*dUzdz;
        Tzz_blk = -2/3.*Visc_blk.*dUxdx - 2/3.*Visc_blk.*dUydy + 4/3.*Visc_blk.*dUzdz;
        Txy_blk =  2.*Visc_blk.*Exy_blk;
        Txz_blk =  2.*Visc_blk.*Exz_blk;
        Tyz_blk =  2.*Visc_blk.*Eyz_blk;
%         % 2D: 
%         Txx_blk        =  4/3.*Visc_blk.*dUxdx - 2/3.*Visc_blk.*dUzdz;
%         Tzz_blk        = -2/3.*Visc_blk.*dUxdx + 4/3.*Visc_blk.*dUzdz;
%         Txz_blk        =  2  .*Visc_blk.*Exz_blk;
        
        %==================================================================
        % STRESS VECTOR, T_i = sigma_ij n_j  (PER UNIT AREA)
        %==================================================================
        Tx = Txx_blk .* xhat_ip_blk + Txy_blk .* yhat_ip_blk + Txz_blk .* zhat_ip_blk;
        Ty = Txy_blk .* xhat_ip_blk + Tyy_blk .* yhat_ip_blk + Tyz_blk .* zhat_ip_blk;
        Tz = Txz_blk .* xhat_ip_blk + Tyz_blk .* yhat_ip_blk + Tzz_blk .* zhat_ip_blk;
        
        %==================================================================
        % NORMAL STRESS, sigma_n = dot(T_i,n_i)
        %==================================================================
        sigma_n = Tx .* xhat_ip_blk + Ty .* yhat_ip_blk + Tz .* zhat_ip_blk;
        
        if strcmp(SETTINGS.is_elastic,'yes')
            %==============================================================
            % VORTICITY
            % (required for stress rotation in NEXT time step)
            %==============================================================
            Vort_x_blk = 0.5.*(sum(Uz_blk.*dNUdy,2) - sum(Uy_blk.*dNUdz,2));
            Vort_y_blk = 0.5.*(sum(Ux_blk.*dNUdz,2) - sum(Uz_blk.*dNUdx,2));
            Vort_z_blk = 0.5.*(sum(Uy_blk.*dNUdx,2) - sum(Ux_blk.*dNUdy,2));
            
%             % 2D:
%             Vort_xz_blk = 0.5.*(  sum(dNUdz.*Vel_blk(:,1:2:end-1),2) ...
%                                 - sum(dNUdx.*Vel_blk(:,2:2:end  ),2) );
            
            % Decay-factor "Xi" for elastic stresses
            ShearG_blk = PHYSICS.ShearG(PhaseID(il:iu));
            Xi         = Visc_blk./(ShearG_blk(:).*dt);
            
            %==============================================================
            % TOTAL STRESSES ARE VISCOUS PLUS OLD ELASTIC STRESSES
            %==============================================================
            Txx_blk    = Txx_blk + Xi.*VAR.Stress_xx(il:iu,ip);
            Tyy_blk    = Tyy_blk + Xi.*VAR.Stress_yy(il:iu,ip);
            Tzz_blk    = Tzz_blk + Xi.*VAR.Stress_zz(il:iu,ip);
            Txy_blk    = Txy_blk + Xi.*VAR.Stress_xy(il:iu,ip);
            Txz_blk    = Txz_blk + Xi.*VAR.Stress_xz(il:iu,ip);
            Tyz_blk    = Tyz_blk + Xi.*VAR.Stress_yz(il:iu,ip);
        end
        
        %==================================================================
        % WRITE DATA INTO GLOBAL STORAGE
        %==================================================================
        VAR.Strain_xx(il:iu,ip) = dUxdx / NUMSCALE.t0;   % strain in 1/s
        VAR.Strain_yy(il:iu,ip) = dUydy / NUMSCALE.t0;   % strain in 1/s
        VAR.Strain_zz(il:iu,ip) = dUzdz / NUMSCALE.t0;   % strain in 1/s
        VAR.Strain_xy(il:iu,ip) = Exy_blk / NUMSCALE.t0; % strain in 1/s
        VAR.Strain_xz(il:iu,ip) = Exz_blk / NUMSCALE.t0; % strain in 1/s
        VAR.Strain_yz(il:iu,ip) = Eyz_blk / NUMSCALE.t0; % strain in 1/s
        
        VAR.Stress_xx(il:iu,ip) = Txx_blk * NUMSCALE.Visc0 / NUMSCALE.t0; % stress in Pa
        VAR.Stress_yy(il:iu,ip) = Tyy_blk * NUMSCALE.Visc0 / NUMSCALE.t0; % stress in Pa
        VAR.Stress_zz(il:iu,ip) = Tzz_blk * NUMSCALE.Visc0 / NUMSCALE.t0; % stress in Pa
        VAR.Stress_xy(il:iu,ip) = Txy_blk * NUMSCALE.Visc0 / NUMSCALE.t0; % stress in Pa
        VAR.Stress_xz(il:iu,ip) = Txz_blk * NUMSCALE.Visc0 / NUMSCALE.t0; % stress in Pa
        VAR.Stress_yz(il:iu,ip) = Tyz_blk * NUMSCALE.Visc0 / NUMSCALE.t0; % stress in Pa
        
        VAR.sigma_n(il:iu,ip)   = sigma_n * NUMSCALE.Visc0 / NUMSCALE.t0; % vertical stress in Pa
        
        VAR.x_ip(il:iu,ip)      = x_ip_blk;
        VAR.y_ip(il:iu,ip)      = y_ip_blk;
        VAR.z_ip(il:iu,ip)      = z_ip_blk;
        
        if strcmp(SETTINGS.is_elastic,'yes')
            VAR.Vort_x(il:iu,ip) = Vort_x_blk;
            VAR.Vort_y(il:iu,ip) = Vort_y_blk;
            VAR.Vort_z(il:iu,ip) = Vort_z_blk;
        end
    end % END OF STRESS EVALUATION LOOP
    
    %======================================================================
    % READJUST START, END AND SIZE OF BLOCK
    %======================================================================
    il = il + nelblk;
    if iblk==nblk-1
        nelblk = nel-iu;
    end
    iu  = iu  + nelblk; % update counter for elements in block
    
end % END OF ELEMENT BLOCK LOOP

end % END FUNCTION calc_stresses_3d