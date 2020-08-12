function [T,WS_SLM] = advect3d_slm_p(MESH,COMM,Vx,Vy,Vz,T,dt,OPTS_SLM,WS_SLM)
% Usage: [T,WS_SLM] = advect3d_slm_p(MESH,COMM,Vx,Vy,Vz,T,dt,OPTS_SLM,WS_SLM)
% 
% Purpose: Performs Semi-Lagrange advection of variable field.
%
% Input:
%   MESH    : [structure] : FE mesh parameters
%   COMM    : [structure] : inter-subdomain communication data
%   WS_SLM  : [structure] : parameters for advection scheme (e.g. els in
%                           which back tracking points have been located
%                           the last time)
%   Vx      : [colvector] : horizontal velocity field
%   Vz      : [colvector] : vertical velocity field
%   T       : [colvector] : variable field to be advected
%   dt      : [scalar]    : time over which is advected
%
% Output:
%   T       : [colvector] : advected variable field
%   WS_SLM  : [structure] : parameters for advection scheme (now els in
%                           which back tracking points have been located)
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
%
                                                                            t0=tic;
% compare_serial_parallel_3d(MESH,COMM,Vx,'Vx',1,0);
% compare_serial_parallel_3d(MESH,COMM,Vy,'Vy',1,0);
% compare_serial_parallel_3d(MESH,COMM,Vz,'Vz',1,0);

% =========================================================================
% Prepare options for backtracking, point search and interpolation
% =========================================================================
OPTS_search.verbose    = 1;
OPTS_search.check_els  = 0;
% OPTS_search.EL2NOD_top = MESH.EL2NOD_top;
% OPTS_search.EL2NOD_bot = MESH.EL2NOD_bot;
OPTS_search.xmin       = MESH.xmin;
OPTS_search.xmax       = MESH.xmax;
OPTS_search.ymin       = MESH.ymin;
OPTS_search.ymax       = MESH.ymax;
OPTS_search.zmin       = MESH.zmin;
OPTS_search.zmax       = MESH.zmax;
OPTS_search.ztol       = 1e-8*MESH.Lz;
OPTS_search.WS         = MESH.WS_tsearch2;
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
    % Make a guess for the elements in which back tracking points are located
    % --> choose one element connected to each node
    [~,ind] = unique(MESH.EL2NOD{1});
    els     = uint32(ceil(ind./size(MESH.EL2NOD{1},1))); clear ind
else
    els     = WS_SLM.els_BT;
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
    case 'RK4' % Runge-Kutta 4th order
        [xyz_BT,els] = backtrack_RK4... % *SUBFUNCTION*
            (MESH,Vx,Vy,Vz,COMM,OPTS_search,OPTS_interp,dt,els);
        iloc                = els>0;
        WS_SLM.els_BT(iloc) = els(iloc);
    case 'PC2' % Predictor-Corrector
        [xyz_BT,els] = backtrack_PC2... % *SUBFUNCTION*
            (MESH,Vx,Vy,Vz,COMM,OPTS_search,OPTS_interp,dt,els);
        iloc                = els>0;
        WS_SLM.els_BT(iloc) = els(iloc);
    case 'EU1' % Euler
        xyz_BT = backtrack_EU1(MESH,Vx,Vy,Vz,OPTS_search,dt); % *SUBFUNCTION*
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
[T_BT,els] = locate_points_interp_vars_3d_p...
    (MESH,COMM,xyz_BT,T,OPTS_search,OPTS_interp,els);
iloc                    = els>0;
WS_SLM.els_INTERP(iloc) = els(iloc);
T(iloc)                 = T_BT(iloc);

% compare_serial_parallel_3d(MESH,COMM,iloc,'iloc',1,0);
% compare_serial_parallel_3d(MESH,COMM,T_BT,'T_BT',1,1);

fprintf(1,'\n SLM ADVECTION (%s,%s)      : %7.2f sec\n',...
    OPTS_SLM.method_BT,OPTS_SLM.method_interp,toc(t0));

end % END OF FUNCTION advect3d_slm_p

% #########################################################################
%                               SUBFUNCTIONS
% #########################################################################

function xyz_BT = backtrack_EU1(MESH,Vx,Vy,Vz,OPTS_search,dt)

% Euler 1st order
% step = [ 1]; step size
% wght = [ 1]; weighting factor
V_BT  = [Vx(:)'; Vy(:)'; Vz(:)'];
xyz_BT = MESH.GCOORD - dt*V_BT;

% Shift backtracking points outsise domain onto boundary
[iout,xyz_BT] = find_pts_outside_domain(MESH.GCOORD,xyz_BT,OPTS_search);
fprintf(' EU1-backtrack: %1i points (located outside domain) shifted onto boundary.\n',...
    length(find(iout)));

end % END OF FUNCTION backtrack_EU1

% #########################################################################

function [xyz_BT,els] = backtrack_PC2(MESH,Vx,Vy,Vz,COMM,OPTS_search,OPTS_interp,dt,els)

% Predictor-Corrector 2nd order
% step = [ 0 .5]; step size
% wght = [ 0  1]; weighting factor

% Gather nodal velocity components
V_nod = [Vx(:) Vy(:) Vz(:)];

% Predictor step
% ==============
xyz_BT = MESH.GCOORD - 0.5*dt*V_nod';

% PTs outside the domain set back to node position for the Corrector step
% =======================================================================
[~,xyz_BT] = find_pts_outside_domain(MESH.GCOORD,xyz_BT,OPTS_search);

% Interpolate velocity at new location
% ====================================
[V_BT,els] = locate_points_interp_vars_3d_p...
    (MESH,COMM,xyz_BT,V_nod,OPTS_search,OPTS_interp,els);

% compare_serial_parallel_3d(MESH,COMM,V_BT(1,:),'Ux_BT',1,0);
% compare_serial_parallel_3d(MESH,COMM,V_BT(2,:),'Uy_BT',1,0);
% compare_serial_parallel_3d(MESH,COMM,V_BT(3,:),'Uz_BT',1,0);

% Corrector step
% ==============
xyz_BT = MESH.GCOORD - dt*V_BT';

% Shift backtracking points outsise domain onto boundary
[iout,xyz_BT] = find_pts_outside_domain(MESH.GCOORD,xyz_BT,OPTS_search);
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