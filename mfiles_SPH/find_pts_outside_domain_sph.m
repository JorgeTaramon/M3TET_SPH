function [iout,thphr_PT] = find_pts_outside_domain_sph(MESH,thphr_PT)
% Usage: [iout,thphr_pt] = find_pts_outside_domain_sph(MESH,thphr_pt)
% 
% Purpose: Finds points that are outside of a shell mesh. If nargout==2
%          shifts these points slightly inside the mesh.
%
% Input:
%   MESH     : [structure] : FE mesh parameters
%   thphr_pt : [matrix]    : spherical coordinates (theta,phi,r) of points
%
% Output:
%   iout     : [vector]    : list of points that were outside the mesh
%   thphr_pt : [matrix]    : new spherical coordinates (theta,phi,r) where
%                            outside points have been radially shifted into
%                            the mesh
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Feb 2016
% JMT Oct 2016: Find outside points using spherical coordinates

verbose   = 1;
nPT       = size(thphr_PT,2);
iout      = true(1,nPT);
th_max    = pi;          % theta max
th_min    = 0;           % theta min
ph_max    = 2*pi;        % phi max
ph_min    = 0;           % phi max
r_cmb     = MESH.r_cmb;  % CMB radius (r min)
r_surf    = MESH.r_surf; % surface radius (r max)
r_tol     = 1e-12 * (r_surf - r_cmb); % Tolerance to place points INSIDE the mesh (rather then exaclty ON the boundary)
th_PT     = thphr_PT(1,:);
ph_PT     = thphr_PT(2,:);
r_PT      = thphr_PT(3,:);
iin       = th_PT >= th_min & th_PT <= th_max  & ...
            ph_PT >= ph_min & ph_PT <= ph_max  & ...
            r_PT  >= r_cmb  + r_tol  & r_PT  <= r_surf - r_tol;
iout(iin) = 0;
ind       = find(iout);
if ~isempty(ind)
    
% % %     % deal with subset of points that have theta < 0 (it means that the
% % %     % points have crossed the North Pole, so, since 0 <= theta <= pi, we
% % %     % need to take the absolute value of theta < 0 and add pi to phi)
% % %     ind_theta_0                  = ind(th_PT(ind)<th_min);
% % %     if ~isempty(ind_theta_0)
% % %         thphr_PT(1,ind_theta_0)  = abs(thphr_PT(1,ind_theta_0));
% % %         thphr_PT(2,ind_theta_0)  = thphr_PT(2,ind_theta_0) + pi;
% % %         % make sure that if we add pi to phi, it does exceed 2pi, otherwise subtract 2pi
% % %         ph_temp                  = thphr_PT(2,ind_theta_0);
% % %         ph_temp(ph_temp >= 2*pi) = ph_temp(ph_temp >= 2*pi) - 2*pi;
% % %         thphr_PT(2,ind_theta_0)  = ph_temp;
% % %     end
% % %     
% % %     % deal with subset of points that have theta > pi (it means that the
% % %     % points have crossed the South Pole, so, since 0 <= theta <= pi, we 
% % %     % need to subtract pi to those theta > pi, then the new theta values 
% % %     % will be pi - (theta > pi) and finally add pi to phi)
% % %     ind_theta_pi                 = ind(th_PT(ind)>th_max);
% % %     if ~isempty(ind_theta_pi)
% % %         delta_th                 = thphr_PT(1,ind_theta_pi) - pi;
% % %         thphr_PT(1,ind_theta_pi) = pi - delta_th;
% % %         thphr_PT(2,ind_theta_pi) = thphr_PT(2,ind_theta_pi) + pi;
% % %         % make sure that if we add pi to phi, it does exceed 2pi, otherwise subtract 2pi
% % %         ph_temp                  = thphr_PT(2,ind_theta_pi);
% % %         ph_temp(ph_temp >= 2*pi) = ph_temp(ph_temp >= 2*pi) - 2*pi;
% % %         thphr_PT(2,ind_theta_pi) = ph_temp;
% % %     end
% % %     
% % %     % deal with subset of points that have phi < 0
% % %     ind_phi_0                    = ind(ph_PT(ind)<ph_min);
% % %     if ~isempty(ind_phi_0)
% % %         thphr_PT(2,ind_phi_0)    = thphr_PT(2,ind_phi_0) + 2*pi;
% % %     end
% % %     
% % %     % deal with subset of points that have phi > 2pi
% % %     ind_phi_2pi                  = ind(ph_PT(ind)>ph_max);
% % %     if ~isempty(ind_phi_2pi)
% % %         thphr_PT(2,ind_phi_2pi)  = thphr_PT(2,ind_phi_2pi) - 2*pi;
% % %     end
    
    % deal with subset of points that are below core-mantle boundary
    ind_cmb                      = ind(r_PT(ind)<=r_cmb + r_tol);
    if ~isempty(ind_cmb)
        thphr_PT(3,ind_cmb)      = r_cmb + r_tol;
    end
    
    % deal with subset of points that are over the surface
    ind_surf                     = ind(r_PT(ind)>=r_surf - r_tol);
    if ~isempty(ind_surf)
        thphr_PT(3,ind_surf)     = r_surf - r_tol;
    end
end

if verbose
    nout = length(find(iout));
    fprintf(' %1i (%.2f%%) points outside of the domain',nout,100*nout/nPT);
    if nargout==2
        fprintf('; shifted onto boundary.\n');
    else
        fprintf('\n');
    end
end

end % END OF FUNCTION find_pts_outside_domain_sph