function [T,WS_SLM] = advect3d_slm_p_sph(MESH,COMM,Ux,Uy,Uz,T,dt,OPTS_SLM,WS_SLM)
% Usage: [T,WS_SLM] = advect3d_slm_p_sph(MESH,COMM,Ux,Uy,Uz,T,dt,OPTS_SLM,WS_SLM)
% 
% Purpose: Performs Semi-Lagrange advection of variable field.
%
% Input:
%   MESH   : [structure] : FE mesh parameters
%   COMM   : [structure] : inter-subdomain communication data
%   WS_SLM : [structure] : parameters for advection scheme (e.g. els in
%                          which back tracking points have been located the
%                          last time)
%   Vth    : [colvector] : theta velocity field
%   Vph    : [colvector] : theta velocity field
%   Vr     : [colvector] : radial velocity field
%   T      : [colvector] : variable field to be advected
%   dt     : [scalar]    : time over which is advected
%
% Output:
%   T      : [colvector] : advected variable field
%   WS_SLM : [structure] : parameters for advection scheme (now els in
%                          which back tracking points have been located)
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Dec 2012
% JH Mar 2013
% JH Aug 2013
% JH Dec 2014
% JMT Oct 2016: Now it works for spherical coordinates
                                                                            t0=tic;
% compare_serial_parallel_3d(MESH,COMM,Vx,'Vx',1,0);
% compare_serial_parallel_3d(MESH,COMM,Vy,'Vy',1,0);
% compare_serial_parallel_3d(MESH,COMM,Vz,'Vz',1,0);

% =========================================================================
% Prepare options for backtracking, point search and interpolation
% =========================================================================
OPTS_search.verbose    = 1;
OPTS_search.check_els  = 0;
OPTS_search.xmin       = MESH.xmin;
OPTS_search.xmax       = MESH.xmax;
OPTS_search.ymin       = MESH.ymin;
OPTS_search.ymax       = MESH.ymax;
OPTS_search.zmin       = MESH.zmin;
OPTS_search.zmax       = MESH.zmax;
OPTS_search.ztol       = 1e-8*MESH.Lz;
% OPTS_search.WS_els_out_cone_no_cross_2pi = MESH.WS_els_out_cone_no_cross_2pi;
% OPTS_search.WS_els_out_cone_cross_2pi    = MESH.WS_els_out_cone_cross_2pi;
% OPTS_search.WS_els_in_cone_no_iso        = MESH.WS_els_in_cone_no_iso;
% OPTS_search.WS_els_in_cone_iso_rot_X_90  = MESH.WS_els_in_cone_iso_rot_X_90;
% OPTS_search.WS_els_in_cone_iso_no_rot    = MESH.WS_els_in_cone_iso_no_rot;
% OPTS_search.WS_els_in_cone_iso_rot_Z_180 = MESH.WS_els_in_cone_iso_rot_Z_180;
OPTS_search.WS_els_in_cone_rot_X_90      = MESH.WS_els_in_cone_rot_X_90;
OPTS_search.WS_els_out_cone              = MESH.WS_els_out_cone;
OPTS_search.WS_els_cross_2pi_rot_Z_180   = MESH.WS_els_cross_2pi_rot_Z_180;
OPTS_search.nthreads   = COMM.nthreads;
% Define interpolation method for velocity at substep coordinates
OPTS_interp.verbose = 1;
switch size(MESH.EL2NOD{1},1)
    case 4
        OPTS_interp.method_interp = 'linear';
    case 10
        OPTS_interp.method_interp = 'quadratic';
    otherwise
        error(' Number of nodes per element must be 4 pr 10');
end

% =========================================================================
% Calculate coordinates of all back tracking points
% =========================================================================
if ~isfield(WS_SLM,'els_BT')
%     % Make a guess for the elements in which back tracking points are located
%     % --> choose one element connected to each node
%     [~,ind] = unique(MESH.EL2NOD{1});
%     els     = uint32(ceil(ind./size(MESH.EL2NOD{1},1))); clear ind
    els = [];
else
    els = WS_SLM.els_BT;
end

% % THIS BLOCK FOR DEBUGGING WORKS !!!
% xyz_BT  = MESH.GCOORD - dt*V_nod;
% EL2NOD  = tetmesh_p2_to_p1(MESH.EL2NOD{1});
% els     = tsearch2(MESH.GCOORD,EL2NOD([1 2 4 3],:),xyz_BT);
% iloc    = els>0;
% lc      = local_coords_3d(MESH.GCOORD,EL2NOD,els(iloc),xyz_BT(:,iloc));
% T(iloc) = interp3d_tet(EL2NOD,els(iloc),lc,T); % interpolate
% return
% % THIS BLOCK FOR DEBUGGING WORKS !!!

switch OPTS_SLM.method_BT
    case 'EU1' % Euler 1st-order
        thphr_BT = backtrack_EU1(MESH,Ux,Uy,Uz,dt); % *SUBFUNCTION*
    case 'PC2' % Predictor-Corrector 2nd-order
        [thphr_BT,els] = backtrack_PC2(MESH,Ux,Uy,Uz,COMM,OPTS_search,OPTS_interp,dt,els); % *SUBFUNCTION*
        iloc                = els>0;
        WS_SLM.els_BT(iloc) = els(iloc);
    case 'RK4' % Runge-Kutta 4th order
        error('RK4 is not working for spherical coordinates (need to be coded)')
        [xyz_BT,els] = backtrack_RK4(MESH,Ux,Uy,Uz,COMM,OPTS_search,OPTS_interp,dt,els); % *SUBFUNCTION*
        iloc                = els>0;
        WS_SLM.els_BT(iloc) = els(iloc);
end

% compare_serial_parallel_3d(MESH,COMM,xyz_BT(1,:)','x_BT',0);
% compare_serial_parallel_3d(MESH,COMM,xyz_BT(2,:)','y_BT',0);
% compare_serial_parallel_3d(MESH,COMM,xyz_BT(3,:)','z_BT',0);

% =========================================================================
% Locate back tracking points in mesh and interpolate temperature
% =========================================================================
if isfield(WS_SLM,'els_INTERP')
    els               = WS_SLM.els_INTERP;
else
    WS_SLM.els_INTERP = els;
end

OPTS_interp.method_interp = OPTS_SLM.method_interp;
OPTS_interp.monotonic     = OPTS_SLM.monotonic;
OPTS_interp.method_wght   = OPTS_SLM.method_wght;
OPTS_interp.verbose       = 1;
[T_BT,els] = locate_points_interp_vars_3d_sph...
    (MESH,COMM,thphr_BT,T,OPTS_search,OPTS_interp,els);
iloc                      = els>0;
WS_SLM.els_INTERP(iloc)   = els(iloc);
T(iloc)                   = T_BT(iloc);

% compare_serial_parallel_3d(MESH,COMM,iloc,'iloc',1,0);
% compare_serial_parallel_3d(MESH,COMM,T_BT,'T_BT',1,1);

fprintf(1,'\n SLM ADVECTION (%s,%s)      : %7.2f sec\n',...
    OPTS_SLM.method_BT,OPTS_SLM.method_interp,toc(t0));

end % END OF FUNCTION advect3d_slm_p_sph

% #########################################################################
%                               SUBFUNCTIONS
% #########################################################################

function thphr_BT = backtrack_EU1(MESH,Vth,Vph,Vr,dt)

% Euler 1st order
% step = [ 1]; step size
% wght = [ 1]; weighting factor

th_BT    = MESH.GCOORD_SPH(1,:) - dt*Vth'./MESH.GCOORD_SPH(3,:); % th_BT = th - dt*v_th_BT/r
ph_BT    = MESH.GCOORD_SPH(2,:) - dt*Vph'./MESH.GCOORD_SPH(3,:); % ph_BT = ph - dt*v_ph_BT/r
r_BT     = MESH.GCOORD_SPH(3,:) - dt*Vr';                        % r_BT  = r  - dt*v_r_BT
thphr_BT = [th_BT; ph_BT; r_BT];

% Shift backtracking points outsise domain onto boundary
[iout,thphr_BT] = find_pts_outside_domain_sph(MESH,thphr_BT);
fprintf(' EU1-backtrack: %1i points (located outside domain) shifted onto boundary.\n',...
    length(find(iout)));

end % END OF FUNCTION backtrack_EU1

% #########################################################################

function [thphr_BT,els] = backtrack_PC2(MESH,Ux,Uy,Uz,COMM,OPTS_search,OPTS_interp,dt,els)
% Predictor-Corrector 2nd order
% step = [ 0 .5]; step size
% wght = [ 0  1]; weighting factor

% theta_cone    = 25*pi/180;     % theta cone in rad
% theta_north   = theta_cone;    % theta angle from the +Z axis to the generatrix
% theta_south   = pi-theta_cone; % theta angle from the -Z axis to the generatrix
% nodes_in_cone = MESH.GCOORD_SPH(1,:) >= theta_south; % boolean vector for nodes in the cone
% % nodes_in_cone = MESH.GCOORD_SPH(1,:) <= theta_north | MESH.GCOORD_SPH(1,:) >= theta_south; % boolean vector for nodes in the cone
% 
% figure(3)
% clf
% plot_coords_cart(MESH.GCOORD_SPH(:,nodes_in_cone),[0 1 1],3)
% figure(4)
% clf
% plot_coords_sph(MESH.GCOORD_SPH(:,nodes_in_cone),[0 1 1],4)

U = [Ux(:)'; Uy(:)'; Uz(:)'];

% Predictor step
% ==============
[thphr_BT,U_perp_unit,r_pred_unit] = predictor_step_with_projections(MESH,U,dt);

% PTs outside the domain set back to node position for the Corrector step
% =======================================================================
[~,thphr_BT] = find_pts_outside_domain_sph(MESH,thphr_BT);

% % % Using the spherical velocities to compute BT points leads to a wrong BT positions (DO NOT USE lines below)  
% % V_nod        = MESH.RR' * U(:); % transform velocity to spherical coordinates
% % V_nod        = reshape(V_nod,3,[]);
% % th_BT        = MESH.GCOORD_SPH(1,:) - 0.5*dt*V_nod(1,:)./MESH.GCOORD_SPH(3,:);                              % th_BT = th - 0.5*dt*v_th/r
% % ph_BT        = MESH.GCOORD_SPH(2,:) - 0.5*dt*V_nod(2,:).*sin(MESH.GCOORD_SPH(1,:))./(MESH.GCOORD_SPH(3,:)); % ph_BT = ph - 0.5*dt*v_ph*sin(th)/r
% % r_BT         = MESH.GCOORD_SPH(3,:) - 0.5*dt*V_nod(3,:);                                                    % r_BT  = r  - 0.5*dt*v_r
% % thphr_BT     = [th_BT; ph_BT; r_BT];
% % [~,thphr_BT] = find_pts_outside_domain_sph(MESH,thphr_BT);

% plot_coords_cart(thphr_BT(:,nodes_in_cone),[1 1 0],3)
% plot_coords_sph(thphr_BT(:,nodes_in_cone),[1 1 0],4)

% Interpolate velocity at new location
% ====================================
[U_BT,els] = locate_points_interp_vars_3d_sph...
    (MESH,COMM,thphr_BT,U',OPTS_search,OPTS_interp,els);

% compare_serial_parallel_3d(MESH,COMM,V_BT(1,:),'Ux_BT',1,0);
% compare_serial_parallel_3d(MESH,COMM,V_BT(2,:),'Uy_BT',1,0);
% compare_serial_parallel_3d(MESH,COMM,V_BT(3,:),'Uz_BT',1,0);

% Corrector step
% ==============
thphr_BT = corrector_step_with_projections(MESH,U_BT',U_perp_unit,r_pred_unit,dt);

% Shift backtracking points outsise domain onto boundary
[iout,thphr_BT] = find_pts_outside_domain_sph(MESH,thphr_BT);

% % % Using the cartesian velocities to compute BT points leads to a wrong BT positions (DO NOT USE lines below) 
% % xyz_BT          = MESH.GCOORD - dt*U_BT';
% % thphr_BT        = cartesian2spherical(xyz_BT);
% % [iout,thphr_BT] = find_pts_outside_domain_sph(MESH,thphr_BT);

% plot_coords_cart(thphr_BT(:,nodes_in_cone),[1 0 0],3)
% plot_coords_sph(thphr_BT(:,nodes_in_cone),[1 0 0],4)
% Ux        = U_BT(:,1)';
% Uy        = U_BT(:,2)';
% Uz        = U_BT(:,3)';
% V_BT      = MESH.RR' * U_BT(:);      % transform velocity to xyz-coordinates
% V_BT      = reshape(V_BT,3,[]);
% Uth       = V_BT(1,:);
% Uph       = V_BT(2,:);
% Ur        = V_BT(3,:);
% plot_velocity_cart(thphr_BT(:,nodes_in_cone),Ux(nodes_in_cone),Uy(nodes_in_cone),Uz(nodes_in_cone),3)
% plot_velocity_sph(thphr_BT(:,nodes_in_cone),Uth(nodes_in_cone),Uph(nodes_in_cone),Ur(nodes_in_cone),4)

fprintf(' PC2-backtrack: %1i points (located outside domain) shifted onto boundary.\n',...
    length(find(iout)));

end % END OF FUNCTION backtrack_PC2

% #########################################################################

function [xyz_BT,els] = backtrack_RK4(MESH,Vx,Vy,Vz,COMM,OPTS_search,OPTS_interp,dt,els)

error(' Verify RK4 first.');

% Runge-Kutta 4th order
step  = [ 0  .5  .5   1 ]; % step size
wght  = [1/6 1/3 1/3 1/6]; % weighting factor
dX_BT = zeros(size(MESH.GCOORD));
xyz_BT = zeros(size(MESH.GCOORD));
V_nod = [Vx(:) Vy(:) Vz(:)];
V_BT  = V_nod; % velocities for back tracking in step "1"
                % are the nodal velocities (since step(1)==0)
iin   = true(1,length(els));
for it=1:4
    % Accumulate backtrack vector
    % ===========================
    dX_BT(:,iin) = dX_BT(:,iin) + wght(it)*dt*V_BT(iin,:)';

    if it<4
        % Calculate new interpolation points
        % ==================================
        xyz_BT(:,iin) = MESH.GCOORD(:,iin) - step(it+1)*dt*V_BT(iin,:)';

        % Check if new interpolation points are outside the domain
        iout_now  = find_pts_outside_domain(MESH.GCOORD,xyz_BT,OPTS_search);
        iin       = iin & ~iout_now;
        
        % PTs outside the domain set back to the node position
        xyz_BT(:,iout_now) = MESH.GCOORD(:,iout_now);
        dX_BT(:,iout_now) = 0;

    else
        % Calculate backtrack points
        % ==========================
        xyz_BT = MESH.GCOORD - dX_BT;
        clear dX_BT
        return
    end

    % Interpolate velocity at new location
    % ====================================
    [V_BT(iin,:),els(iin)] = locate_points_interp_vars_3d_p...
        (MESH,COMM,xyz_BT(:,iin),V_nod,OPTS_search,OPTS_interp,els(iin));
end

% Shift backtracking points outsise domain onto boundary
[iout,xyz_BT] = find_pts_outside_domain(MESH.GCOORD,xyz_BT,OPTS_search);
fprintf(' RK4-backtrack: %1i points (located outside domain) shifted onto boundary.\n',...
    length(find(iout)));

end % END OF FUNCTION backtrack_RK4

function [thphr_BT,U_perp_unit,r_BT_unit] = predictor_step_with_projections(MESH,U,dt)
% Usage: [thphr_BT,U_perp_unit,r_BT_unit] = predictor_step_with_projections(MESH,U,dt)
%
% Purpose: Compute the backtrack points for predictor step. Parallel and
%          perpendicular projections for the velocity are computed.
%
% Input:
%   MESH        : [structure] : FE mesh
%   U           : [matrix]    : velocity field in Cartesian coordinates
%   dt          : [scalar]    : time over which is advected
%
% Output:
%   thphr_BT    : [matrix] : spherical coordinates for BT points
%   U_perp_unit : [matrix] : unit U_perpendicular vector
%   r_BT_unit   : [matrix] : unit radial vector for BT points
%                                                                                  X r_BT_unit
%                                                                                 /
%                                                                                /
%                                                                               / 
%                                                                              /  
%                                                                           _ X r'
%                      X r_unit                                        __--  /
%                       \                                        __ --      /
%                        \                                 __ --           /
%                         \                       l' __ --                /
%                          \  U_par            __ --                     /
%                           X             _ -- _______                  /
%               U            \        ___ ----         ---- ___        /
%                 X -- ___    \  _ --                           -- _  /
%                          --- X   l_perp = |U_perp|*step*dt = r*da  X <-- = r_BT_unit * |r|
%                      __ --  r \                                   /
%                   X-       /   \                                 /
%                 U_perp     \    \                               /
%                             \    \                             /
%                              \    \                           /
%       l_par = |U_par|*step*dt \    \                         /
%                                \    \                       /
%                                 ---  +                     X  r_BT = r_BT_unit * |r - l_par|
%                                       \                   /     \\
%                                        \                 /       \\
%                                         \               /         \\
%                                          \             /       ====================================
%                                           \           /       || THIS IS THE POINT TO BE COMPUTED || 
%                                            \    da   /         ====================================
%             dot(U,r_unit)                   \  ___  /
%   U_par  = -------------- r_unit             \/   \/ 
%              |r_unit|^2                       \   /
%                                                \ /
%   U_perp = U - U_par                            o
%
%   1st STEP: Move the point in the perpendicular direction to the unit vector
%       
%       In order to get the unit vector for the BT point (r_BT_unit), geometrically: 
%
%                           l'
%                           - = tg(da)
%                           r
%                                                  l_perp
%       On the other hand: l_perp = r*da --> da = -------- , so:
%                                                    r
%                            ( l_perp )
%                   l' = r*tg|--------|
%                            (   r    )
%
%       Then, the r' vector is given by:
%
%                              U_perp
%                   r' = r - ---------- * l'
%                             |U_perp|
%
%       and its unit vector is r_BT_unit
%
%   2nd STEP: Move the point in the parallel direction to the unit vector (of the BT point --> r_BT_unit) 
%       
%       The coordinates of the BT point are given by:
%       
%                   r_BT = r_BT_unit * |r - l_par|
%

r_unit         = MESH.GCOORD./repmat(MESH.GCOORD_SPH(3,:),3,1);                                   % unit radial vector for the points
U_par          = (repmat(dot(U,r_unit),3,1) .* r_unit)./repmat(sqrt(sum(abs(r_unit).^2)).^2,3,1); % parallel velocity component (vector) to the unit vector
U_perp         = U - U_par;                                % perpendicular velocity component (vector) to the unit vector
U_perp_unit    = U_perp./repmat(sqrt(sum(U_perp.^2)),3,1); % unit vector for the perpendicular direction
if sum(any(isnan(U_perp_unit))) > 0
    % This may happen because:
    % 1) Some nodes have velocity = 0 (they are fixed). This might happen, for example when we use Free Slip BC and we need to fix one or two nodes to avoid net rotation.
    %    The velocity U = [0; 0; 0] and therefore U_par and U_perp are [0; 0; 0] 
    % 2) The unit vector (r_unit) and the velocity vector (U) are in the same direction, meaning that there is only parallel component (U_par = U and U_perp = [0; 0; 0]).
    % Thus,for both cases, when computing U_perp_unit = U_perp/|U_perp| it will give NaNs since it is U_perp_unit would be 0/0 !!
    ind_NaNs = find(sum(isnan(U_perp_unit)));
    % To avoid NaNs:
    for i = 1:size(ind_NaNs,2)
        if sum(U(:,ind_NaNs(i))) == 0 || isequal(U(:,ind_NaNs(i)),U_par(:,ind_NaNs(i)))
            % Check if 
            %   1. NaNs positions in U_perp_unit correspond to the same positions of zero velocity in U or
            %   2. The parallel velocity component (U_par) is the same than U vector 
            % If so, change NaNs values in U_perp_unit by zeros.
            U_perp_unit(:,ind_NaNs(i)) = zeros(3,1);
        else
            error(' This should not happen');
        end
    end
end
l_par          = sqrt(sum(U_par.^2)) * 0.5 * dt;           % parallel length (scalar) to the unit vector for the BT points
l_perp         = sqrt(sum(U_perp.^2)) * 0.5 * dt;          % arc-length (scalar) in the perpendicular direction to the unit vector for the BT points
l_prime        = MESH.GCOORD_SPH(3,:).* tan(l_perp./MESH.GCOORD_SPH(3,:)); % straight length (l') from the point r to r'
l_prime_vector = U_perp_unit.*repmat(l_prime,3,1);                         % l' projected on the U_perp direction
r_prime        = MESH.GCOORD - l_prime_vector;                             % coordinates of r'
r_BT_unit      = r_prime./repmat(sqrt(sum(abs(r_prime).^2)),3,1);          % unit radial vector for the BT points
xyz_BT         = r_BT_unit.*repmat(abs(MESH.GCOORD_SPH(3,:) - l_par),3,1); % Cartesian coordinates for the BT points
thphr_BT       = cartesian2spherical(xyz_BT);                              % spherical coordinates for the BT points

end % END OF SUBFUNCTION predictor_step_with_projections

function thphr_BT = corrector_step_with_projections(MESH,U,U_perp_unit,r_pred_unit,dt)
% Usage: thphr_BT = corrector_step_with_projections(MESH,U,U_perp_unit,r_pred_unit,dt)
%
% Purpose: Compute the backtrack points for corrector step. Parallel and
%          perpendicular projections for the velocity are computed.
%
% Input:
%   MESH        : [structure] : FE mesh
%   U           : [matrix]    : velocity field in Cartesian coordinates
%   U_perp_unit : [matrix]    : unit U_perpendicular vector (in the
%                               original points)
%   r_pred_unit : [matrix]    : unit radial vector for predictor points
%   dt          : [scalar]    : time over which is advected
%
% Output:
%   thphr_BT    : [matrix] : spherical coordinates for BT points

U_par          = (repmat(dot(U,r_pred_unit),3,1) .* r_pred_unit)./repmat(sqrt(sum(abs(r_pred_unit).^2)).^2,3,1); % parallel velocity component (vector) to the unit vector in the predictor step
U_perp         = U - U_par;                 % perpendicular velocity component (vector) to the unit vector in the predictor step
l_par          = sqrt(sum(U_par.^2)) * dt;  % parallel length (scalar) to the unit vector for the BT points
l_perp         = sqrt(sum(U_perp.^2)) * dt; % arc-length (scalar) in the perpendicular direction to the unit vector for the BT points
l_prime        = MESH.GCOORD_SPH(3,:).* tan(l_perp./MESH.GCOORD_SPH(3,:)); % straight length (l') from the point r to r'
l_prime_vector = U_perp_unit.*repmat(l_prime,3,1);                         % l' projected on the U_perp direction of the original points
r_prime        = MESH.GCOORD - l_prime_vector;                             % coordinates of r'
r_BT_unit      = r_prime./repmat(sqrt(sum(abs(r_prime).^2)),3,1);          % unit vector for the BT points
xyz_BT         = r_BT_unit.*repmat(abs(MESH.GCOORD_SPH(3,:) - l_par),3,1); % Cartesian coordinates for the BT points
thphr_BT       = cartesian2spherical(xyz_BT);                              % spherical coordinates for the BT points

end % END OF SUBFUNCTION corrector_step_with_projections

function plot_coords_cart(gTH_PT,color_nodes,Fig_num)

figure(Fig_num)
gX_PT = spherical2cartesian(gTH_PT);
% Plot the shell
lightGrey           = 0.90*[1 1 1]; % colour for the shell boundaries
[x_sph,y_sph,z_sph] = sphere(20);
x_sph               = x_sph*3471;
y_sph               = y_sph*3471;
z_sph               = z_sph*3471;
axis equal
surface(x_sph,y_sph,z_sph,'FaceColor','none','EdgeColor',lightGrey)
[x_sph,y_sph,z_sph] = sphere(30);
x_sph               = x_sph*6371;
y_sph               = y_sph*6371;
z_sph               = z_sph*6371;
surface(x_sph,y_sph,z_sph,'FaceColor','none','EdgeColor',lightGrey)
axis([-6371 6371 -6371 6371 -6371 6371])
view(142.5,30)
grid on
hold on

scatter3(gX_PT(1,:),gX_PT(2,:),gX_PT(3,:),20,...
    'MarkerEdgeColor','k','MarkerFaceColor',color_nodes)
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
drawnow

end % END OF SUBFUNCTION plot_coords_cart

function plot_coords_sph(gTH_PT,color_nodes,Fig_num)

figure(Fig_num)
scatter3(gTH_PT(1,:),gTH_PT(2,:),gTH_PT(3,:)/1000,20,'MarkerEdgeColor','k','MarkerFaceColor',color_nodes)
axis equal
axis([0 pi 0 2*pi 3 6.5])
view(142.5,30)
grid on
hold on
xlabel('\theta (rad)')
ylabel('\phi (rad)')
zlabel('r (10^3 km)')
drawnow

end % END OF SUBFUNCTION plot_coords_sph

function plot_velocity_cart(gTH_PT,Ux,Uy,Uz,Fig_num)

figure(Fig_num)
gX_PT = spherical2cartesian(gTH_PT);
% Plot the shell
lightGrey           = 0.90*[1 1 1]; % colour for the shell boundaries
[x_sph,y_sph,z_sph] = sphere(20);
x_sph               = x_sph*3471;
y_sph               = y_sph*3471;
z_sph               = z_sph*3471;
axis equal
surface(x_sph,y_sph,z_sph,'FaceColor','none','EdgeColor',lightGrey)
[x_sph,y_sph,z_sph] = sphere(30);
x_sph               = x_sph*6371;
y_sph               = y_sph*6371;
z_sph               = z_sph*6371;
surface(x_sph,y_sph,z_sph,'FaceColor','none','EdgeColor',lightGrey)
axis([-6371 6371 -6371 6371 -6371 6371])
view(142.5,30)
grid on
hold on
quiver3(gX_PT(1,:),gX_PT(2,:),gX_PT(3,:),Ux,Uy,Uz,2,'k')
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
drawnow

end % END OF SUBFUNCTION plot_velocity_cart

function plot_velocity_sph(gTH_PT,Uth,Uph,Ur,Fig_num)

figure(Fig_num)
quiver3(gTH_PT(1,:),gTH_PT(2,:),gTH_PT(3,:)/1000,Uth./gTH_PT(3,:),Uph./gTH_PT(3,:),Ur,1,'k')
axis equal
axis([0 pi 0 2*pi 3 6.5])
view(142.5,30)
grid on
hold on
xlabel('\theta (rad)')
ylabel('\phi (rad)')
zlabel('r (10^3 km)')
drawnow

end % END OF SUBFUNCTION plot_velocity_sph