function [lc,delS] = compute_lc_surface_curved(GCOORD_SPH_surf,EL2NOD_face,curved_els,gTH_PT)
% Usage: [lc,delS] = compute_lc_surface_curved(GCOORD_SPH_surf,EL2NOD_face,curved_els,gTH_PT)
%
% Purpose:
%   Returns local coordinates of points "gTH_PT" in curved edge elements
%   "els". 
%
% Input:
%   GCOORD_SPH_surf : [matrix]    : colatitude and longitude (theta,phi) 
%                                   (in rads) of surface nodes of the mesh
%   EL2NOD_face     : [matrix]    : connectivity matrix for faces on the
%                                   surface
%   curved_els      : [rowvector] : curved edge elements in which each
%                                   point is located 
%   gTH_PT          : [matrix]    : coordinates of points to be located
%
% Output:
%   lc              : [matrix]    : local coordinates of points in each
%                                   curved edge element
%   delS            : [matrix]    : 
%
% Part of 3D convection code M3TET_MG, which is a simplified version of the
% parallel 3D mantle convection code M3TET (developed by
% J.Hasenclever & J.Phipps Morgan, 2007-2010)
% Email contact: joerg.hasenclever@zmaw.de
% For numerical methods see online Ph.D. thesis
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)

% JH March 2011, JPM June 2011
% JMT Oct 2016 : works for spherical coordinates

FigNo = 0;

if size(GCOORD_SPH_surf,1) == 3
    GCOORD_SPH_surf(3,:) = [];
end
if size(gTH_PT,1) == 3
    gTH_PT(3,:) = [];
end

method = 'std';
% method = 'opt';

ns = length(curved_els); % Number of points
th  = reshape(GCOORD_SPH_surf(1,EL2NOD_face(1:3,curved_els)),3,ns); % Coords of Vertex nodes
ph  = reshape(GCOORD_SPH_surf(2,EL2NOD_face(1:3,curved_els)),3,ns);

thp = gTH_PT(1,:); % Global coords of points in each element
php = gTH_PT(2,:);

% Definition of local coordinates
% Make sure that shape (interpolation) functions are defined accordingly!
%
%   / \    3
%    |     | \
% lc(2,:)  |   \
%    |     |     \
%    |     1 - - - 2     
%        -- lc(1,:) -->
lc      = zeros(2,length(curved_els)); % use linear estimates as start guesses
denom   = -th(1,:).*ph(3,:)+th(1,:).*ph(2,:)-th(2,:).*ph(1,:)+th(2,:).*ph(3,:)+th(3,:).*ph(1,:)-th(3,:).*ph(2,:);
lc(1,:) = -(-th(1,:).*php+th(1,:).*ph(3,:)-th(3,:).*ph(1,:)+thp.*ph(1,:)+th(3,:).*php-thp.*ph(3,:)) ./ denom;
lc(2,:) =  (th(1,:).*ph(2,:)-th(1,:).*php-th(2,:).*ph(1,:)+th(2,:).*php+thp.*ph(1,:)-thp.*ph(2,:)) ./ denom;

delTH_tol = 1e-8;
nimprove  = 6; % iterative improvement
for improve = 1:nimprove
    % get local estimate of dTH/dS matrix at each guess location
    % (1) get local shapefunctions and derivatives
    [N6,dN6] = shapefun_6nodeVEC(lc(1,:),lc(2,:)); % SUBFUNCTION
    [N6_v2,dN6_v2] = sf_dsf_tri36712([lc(1,:);lc(2,:)],6,'cell','std');
    
    % find misfit delTH between current real position and guess position
    th0 = sum( reshape(GCOORD_SPH_surf(1,EL2NOD_face(:,curved_els))',6,[]).*N6, 1);
    ph0 = sum( reshape(GCOORD_SPH_surf(2,EL2NOD_face(:,curved_els))',6,[]).*N6, 1);
    TH0 = [th0; ph0];

    delTH = gTH_PT - TH0;
    err   = sqrt(sum(delTH .* delTH,1));
    
    [max_err,imax_err] = max(err);
    
    figure(994);clf
    plot(gTH_PT(1,imax_err),gTH_PT(2,imax_err),'xr')
    axis equal
    hold on
    trimesh(EL2NOD_face(:,curved_els(imax_err))',GCOORD_SPH_surf(1,:)',GCOORD_SPH_surf(2,:)','Color',[0 0 0])
    plot(TH0(1,imax_err),TH0(2,imax_err),'ob')
    plot(GCOORD_SPH_surf(1,EL2NOD_face(1,curved_els(imax_err))),GCOORD_SPH_surf(2,EL2NOD_face(1,curved_els(imax_err))),'^r')
    plot(GCOORD_SPH_surf(1,EL2NOD_face(2,curved_els(imax_err))),GCOORD_SPH_surf(2,EL2NOD_face(2,curved_els(imax_err))),'or')
    plot(GCOORD_SPH_surf(1,EL2NOD_face(3,curved_els(imax_err))),GCOORD_SPH_surf(2,EL2NOD_face(3,curved_els(imax_err))),'dr')
    plot(GCOORD_SPH_surf(1,EL2NOD_face(4,curved_els(imax_err))),GCOORD_SPH_surf(2,EL2NOD_face(4,curved_els(imax_err))),'^k')
    plot(GCOORD_SPH_surf(1,EL2NOD_face(5,curved_els(imax_err))),GCOORD_SPH_surf(2,EL2NOD_face(5,curved_els(imax_err))),'ok')
    plot(GCOORD_SPH_surf(1,EL2NOD_face(6,curved_els(imax_err))),GCOORD_SPH_surf(2,EL2NOD_face(6,curved_els(imax_err))),'dk')
    
    if max_err<delTH_tol
        break
    end
    
    if FigNo
%         [~,ipt_plot] = sort(err,'descend');
%         ipt_plot     = ipt_plot(1:10);
        if improve==1
            ipt_plot     = imax_err;
        end
        if improve == 1
            figure(FigNo);clf
        end
        subplot(3,1,1); plot(improve,log10(err(ipt_plot)),'ro'); hold on
        subplot(3,1,2); plot(improve,th0(ipt_plot),'ro'); hold on
        subplot(3,1,3); plot(improve,ph0(ipt_plot),'ro'); hold on
    end
    
    if strcmp(method,'std')
        delS = zeros(size(delTH));
        for iel=1:length(th0)
            dN6_iel     = dN6(:,:,iel);
            delTH_iel   = delTH(:,iel);
            gTH_iel     = GCOORD_SPH_surf(:,EL2NOD_face(:,curved_els(iel)));
            dTHdS       = gTH_iel * dN6_iel'; % outer product
            delS_iel    = dTHdS \ delTH_iel;
            delS(:,iel) = delS_iel;
        end
    else % then use  strcmp(method,'opt')
        % (2) make matrix terms dTH/dS = [sum(GCOORD_SPH .* dN6)]
        dTHdS2     = zeros(2,2,length(th0)); % ns 2x2 matrices, one for each point
        inv_dTHdS = zeros(2,2,length(th0)); % ns 2x2 matrices, one for each point
        % dth/dS terms
        dTHdS2(1,1,:) = sum( reshape(GCOORD_SPH_surf(1,EL2NOD_face(:,curved_els))',6,[]).*reshape(dN6(1,:,:),6,[]), 1);
        dTHdS2(1,2,:) = sum( reshape(GCOORD_SPH_surf(1,EL2NOD_face(:,curved_els))',6,[]).*reshape(dN6(2,:,:),6,[]), 1);
        % dph/dS terms    
        dTHdS2(2,1,:) = sum( reshape(GCOORD_SPH_surf(2,EL2NOD_face(:,curved_els))',6,[]).*reshape(dN6(1,:,:),6,[]), 1);
        dTHdS2(2,2,:) = sum( reshape(GCOORD_SPH_surf(2,EL2NOD_face(:,curved_els))',6,[]).*reshape(dN6(2,:,:),6,[]), 1);

        % determinants of the dTHdS matrices
        det_dTHdS = dTHdS2(1,1,:).*dTHdS2(2,2,:) ...
                  - dTHdS2(1,2,:).*dTHdS2(2,1,:);
        
        % 1/det(dXdSs)
        invdet = 1./det_dTHdS;
        
        inv_dTHdS(1,1,:) = invdet .* dTHdS2(2,2,:);
        inv_dTHdS(2,1,:) = invdet .* dTHdS2(1,2,:);
        
        inv_dTHdS(1,2,:) = invdet .* dTHdS2(2,1,:);
        inv_dTHdS(2,2,:) = invdet .* dTHdS2(1,1,:);

        delS2(1,:) =  inv_dTHdS(1,1,:) .* reshape(delTH(1,:),1,1,[]) ...
                   + inv_dTHdS(1,2,:) .* reshape(delTH(2,:),1,1,[]);
        delS2(2,:) =  inv_dTHdS(2,1,:) .* reshape(delTH(1,:),1,1,[]) ...
                   + inv_dTHdS(2,2,:) .* reshape(delTH(2,:),1,1,[]);
    end

    lc = lc + delS; % update guess
end

end % END OF FUNCTION local_coords_curved_3d_sph

% ##############################################################################################################################################################################################
%                               SUBFUNCTIONS
% ##############################################################################################################################################################################################

function [N6,dN6] = shapefun_6nodeVEC(r,s)
% Find shape functions and their derivatives at given points on the
% master element for a 6 node triangle

%        3
%        | \
% s-axis 6   5
%        |     \
%        1 - 4 - 2
%          r-axis

% Node  r  s
%  1    0  0
%  2    1  0
%  3    0  1
%  4   .5  0
%  5   .5 .5
%  6    0 .5

npts = length(r);

N6  = zeros(6,npts);
dN6 = zeros(2,6,npts);
t    = 1-r-s;

% 3 vertex nodes
N6( 1,:) = t.*(2.*t-1);
N6( 2,:) = r.*(2.*r-1);
N6( 3,:) = s.*(2.*s-1);
% 3 edge nodes
N6( 4,:) = 4.*r.*t;
N6( 5,:) = 4.*r.*s;
N6( 6,:) = 4.*s.*t;

% derivatives (3 for each node)
dN6(1,1,:) = -4*t+1;   % dN1/dr
dN6(2,1,:) = -4*t+1;   % dN1/ds

dN6(1,2,:) =  4*r-1;   % dN2/dr
dN6(2,2,:) =  0;       % dN2/ds

dN6(1,3,:) =  0;       % dN3/dr
dN6(2,3,:) =  4*s-1;   % dN3/ds

dN6(1,4,:) =  4*(t-r); % dN4/dr
dN6(2,4,:) = -4*r;     % dN4/ds

dN6(1,5,:) =  4*s;     % dN5/dr
dN6(2,5,:) =  4*r;     % dN5/ds

dN6(1,6,:) = -4*s;     % dN6/dr
dN6(2,6,:) =  4*(t-s); % dN6/ds

end % END OF FUNCTION shapefun_6nodeVEC