function dt = calculate_dt(MESH,COMM,SETTINGS,NUMSCALE,VAR)
% Usage: dt = calculate_dt(MESH,COMM,SETTINGS,NUMSCALE,VAR)
% 
% Purpose: Calculate next time step based on (1) element size and thermal
%          diffusivity and (2) maximum advection distance defined.
%
% Input:
%   MESH     : [structure] : FE mesh parameters
%   COMM     : [structure] : communication data
%   SETTINGS : [structure] : model parameters
%   NUMSCALE : [structure] : numerical scaling parameters
%   PHYSICS  : [structure] : physical properties
%   VAR      : [structure] : major variable fields
%
% Output:
%   dt       : [scalar]    : length of next time step
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Dec 2012
% JH March 2013
% JH Dec 2014 : simplified
%

% % "Length" (node distance) for each element
% ellength = MESH.len_el;
% 
% % Diffusion limit for time step (Courant criterium)
% dt_diff = min(ellength).^2 / (2*PHYSICS.kappa);
% if isfield(SETTINGS,'dt_diff_limit')
%     dt_diff = dt_diff * SETTINGS.dt_diff_limit;
% end

% Advection limit for time step
maxUx = max(abs(VAR.Ux));
maxUy = max(abs(VAR.Uy));
maxUz = max(abs(VAR.Uz));
if maxUx==0 && maxUy==0 && maxUz==0
    dt_adv   = 0;
    dt       = SETTINGS.dtmax;
else
    dt_adv_x = MESH.Lx/maxUx;
    dt_adv_y = MESH.Ly/maxUy;
    dt_adv_z = MESH.Lz/maxUz;
    dt_adv   = SETTINGS.dt_adv_limit * min([dt_adv_x dt_adv_y dt_adv_z]);
    dt       = min(SETTINGS.dtmax,dt_adv);
end

if COMM.nsd>1
    dt = COMM.minLabs(dt);
end

fprintf(SETTINGS.fid_log,...
        ' TIME STEP LIMITS [%s]: dt(adv)=%0.1e  new dt=%0.1e\n',...
        NUMSCALE.unit_t,dt_adv,dt);

end % END OF FUNCTION calculate_dt