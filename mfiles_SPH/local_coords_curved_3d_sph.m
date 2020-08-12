function [lc,delS] = local_coords_curved_3d_sph(GCOORD_SPH,EL2NOD,curved_els,gTH_PT)
% Usage: [lc,delS] = local_coords_curved_3d_sph(GCOORD_SPH,EL2NOD,curved_els,gTH_PT)
%
% Purpose:
%   Returns local coordinates of points "gTH_PT" in curved edge elements
%   "els". 
%
% Input:
%   GCOORD_SPH : [matrix]    : spherical coodinates (theta,phi,r) 
%                              (in rad, rad and km) of the mesh
%   EL2NOD     : [matrix]    : connectivity matrix for GCOORD_SPH
%   curved_els : [rowvector] : curved edge elements in which each point is
%                              located 
%   gTH_PT     : [matrix]    : coordinates of points to be located
%
% Output:
%   lc         : [matrix]    : local coordinates of points in each curved
%                              edge element
%   delS       : [matrix]    : 
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

% method = 'std';
method = 'opt';

ns = length(curved_els); % Number of points
th  = reshape(GCOORD_SPH(1,EL2NOD(1:4,curved_els)),4,ns); % Coords of Vertex nodes
ph  = reshape(GCOORD_SPH(2,EL2NOD(1:4,curved_els)),4,ns);
r   = reshape(GCOORD_SPH(3,EL2NOD(1:4,curved_els)),4,ns);

thp = gTH_PT(1,:); % Global coords of points in each element
php = gTH_PT(2,:);
rp  = gTH_PT(3,:);

lc = zeros(3,length(curved_els)); % use linear estimates as start guesses

lc(1,:) =   (  (th(2,:)-thp)    .*(ph(3,:)-php)    .*(r(4,:) -rp)      + (th(3,:)-thp)    .*(ph(4,:)-php)    .*(r(2,:) -rp)      + (th(4,:)-thp)    .*(ph(2,:)-php)    .*(r(3,:) -rp)        ...
             - (r(2,:) -rp)     .*(ph(3,:)-php)    .*(th(4,:)-thp)     - (r(3,:) -rp)     .*(ph(4,:)-php)    .*(th(2,:)-thp)     - (r(4,:) -rp)     .*(ph(2,:)-php)    .*(th(3,:)-thp) )     ...
         ./ (  (th(2,:)-th(1,:)).*(ph(3,:)-ph(1,:)).*(r(4,:) -r(1,:))  + (th(3,:)-th(1,:)).*(ph(4,:)-ph(1,:)).*(r(2,:) -r(1,:))  + (th(4,:)-th(1,:)).*(ph(2,:)-ph(1,:)).*(r(3,:) -r(1,:))    ...
             - (r(2,:) -r(1,:)) .*(ph(3,:)-ph(1,:)).*(th(4,:)-th(1,:)) - (r(3,:) -r(1,:)) .*(ph(4,:)-ph(1,:)).*(th(2,:)-th(1,:)) - (r(4,:) -r(1,:)) .*(ph(2,:)-ph(1,:)).*(th(3,:)-th(1,:)) );

lc(2,:) =   (  (th(1,:)-thp)    .*(ph(3,:)-php)    .*(r(4,:) -rp)      + (th(3,:)-thp)    .*(ph(4,:)-php)    .*(r(1,:) -rp)      + (th(4,:)-thp)    .*(ph(1,:)-php)    .*(r(3,:) -rp)        ...
             - (r(1,:)-rp)      .*(ph(3,:)-php)    .*(th(4,:)-thp)     - (r(3,:) -rp)     .*(ph(4,:)-php)    .*(th(1,:)-thp)     - (r(4,:) -rp)     .*(ph(1,:)-php)    .*(th(3,:)-thp) )     ...
         ./ (  (th(1,:)-th(2,:)).*(ph(3,:)-ph(2,:)).*(r(4,:) -r(2,:))  + (th(3,:)-th(2,:)).*(ph(4,:)-ph(2,:)).*(r(1,:) -r(2,:))  + (th(4,:)-th(2,:)).*(ph(1,:)-ph(2,:)).*(r(3,:) -r(2,:))    ...
             - (r(1,:)-r(2,:))  .*(ph(3,:)-ph(2,:)).*(th(4,:)-th(2,:)) - (r(3,:) -r(2,:)) .*(ph(4,:)-ph(2,:)).*(th(1,:)-th(2,:)) - (r(4,:) -r(2,:)) .*(ph(1,:)-ph(2,:)).*(th(3,:)-th(2,:)) );
             
lc(3,:) =   (  (th(1,:)-thp)    .*(ph(2,:)-php)    .*(r(4,:) -rp)      + (th(2,:)-thp)    .*(ph(4,:)-php)    .*(r(1,:) -rp)      + (th(4,:)-thp)    .*(ph(1,:)-php)    .*(r(2,:) -rp)        ...
             - (r(1,:)-rp)      .*(ph(2,:)-php)    .*(th(4,:)-thp)     - (r(2,:) -rp)     .*(ph(4,:)-php)    .*(th(1,:)-thp)     - (r(4,:) -rp)     .*(ph(1,:)-php)    .*(th(2,:)-thp) )     ...
         ./ (  (th(1,:)-th(3,:)).*(ph(2,:)-ph(3,:)).*(r(4,:) -r(3,:))  + (th(2,:)-th(3,:)).*(ph(4,:)-ph(3,:)).*(r(1,:) -r(3,:))  + (th(4,:)-th(3,:)).*(ph(1,:)-ph(3,:)).*(r(2,:) -r(3,:))    ...
             - (r(1,:)-r(3,:))  .*(ph(2,:)-ph(3,:)).*(th(4,:)-th(3,:)) - (r(2,:) -r(3,:)) .*(ph(4,:)-ph(3,:)).*(th(1,:)-th(3,:)) - (r(4,:) -r(3,:)) .*(ph(1,:)-ph(3,:)).*(th(2,:)-th(3,:)) );

delTH_tol = 1e-8;
nimprove  = 6; % iterative improvement
for improve = 1:nimprove
    % get local estimate of dTH/dS matrix at each guess location
    % (1) get local shapefunctions and derivatives
    [N10,dN10] = shapefun_10nodeVEC(lc(1,:),lc(2,:),lc(3,:)); % SUBFUNCTION

    % find misfit delTH between current real position and guess position
    th0 = sum( reshape(GCOORD_SPH(1,EL2NOD(:,curved_els))',10,[]).*N10, 1);
    ph0 = sum( reshape(GCOORD_SPH(2,EL2NOD(:,curved_els))',10,[]).*N10, 1);
    r0  = sum( reshape(GCOORD_SPH(3,EL2NOD(:,curved_els))',10,[]).*N10, 1);
    TH0 = [th0;...
           ph0;...
           r0];

    delTH = gTH_PT - TH0;
    err   = sqrt(sum(delTH .* delTH,1));
    
    [max_err,imax_err] = max(err);
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
        subplot(4,1,1); plot(improve,log10(err(ipt_plot)),'ro'); hold on
        subplot(4,1,2); plot(improve,th0(ipt_plot),'ro'); hold on
        subplot(4,1,3); plot(improve,ph0(ipt_plot),'ro'); hold on
        subplot(4,1,4); plot(improve,r0(ipt_plot),'ro'); hold on
    end
    
    if strcmp(method,'std')
        delS = zeros(size(delTH));
        for iel=1:length(th0)
            dN10_iel  = dN10(:,:,iel);
            delTH_iel = delTH(:,iel);
            gTH_iel   = GCOORD_SPH(:,EL2NOD(:,curved_els(iel)));
    %         % diagonal terms of dXdS (for checking)
    %         thdN10    = gTH_iel(1,:) * dN10_iel(1,:)';
    %         phdN10    = gTH_iel(2,:) * dN10_iel(2,:)';
    %         rdN10     = gTH_iel(3,:) * dN10_iel(3,:)';
            dTHdS     = gTH_iel * dN10_iel'; % outer product
%             J         = dN10_iel*gTH_iel'; % in assembly
            delS_iel = dTHdS \ delTH_iel;
            delS(:,iel) = delS_iel;
        end
    else % then use  strcmp(method,'opt')
        % (2) make matrix terms dTH/dS = [sum(GCOORD_SPH .* dN10)]
        % dth/dS terms
        dTHdS     = zeros(3,3,length(th0)); % ns 3x3 matrices, one for each point
        inv_dTHdS = zeros(3,3,length(th0)); % ns 3x3 matrices, one for each point

        dTHdS(1,1,:) = sum( reshape(GCOORD_SPH(1,EL2NOD(:,curved_els))',10,[]).*reshape(dN10(1,:,:),10,[]), 1);
        dTHdS(1,2,:) = sum( reshape(GCOORD_SPH(1,EL2NOD(:,curved_els))',10,[]).*reshape(dN10(2,:,:),10,[]), 1);
        dTHdS(1,3,:) = sum( reshape(GCOORD_SPH(1,EL2NOD(:,curved_els))',10,[]).*reshape(dN10(3,:,:),10,[]), 1);
        % dph/dS terms    
        dTHdS(2,1,:) = sum( reshape(GCOORD_SPH(2,EL2NOD(:,curved_els))',10,[]).*reshape(dN10(1,:,:),10,[]), 1);
        dTHdS(2,2,:) = sum( reshape(GCOORD_SPH(2,EL2NOD(:,curved_els))',10,[]).*reshape(dN10(2,:,:),10,[]), 1);
        dTHdS(2,3,:) = sum( reshape(GCOORD_SPH(2,EL2NOD(:,curved_els))',10,[]).*reshape(dN10(3,:,:),10,[]), 1);
        % dr/dS terms
        dTHdS(3,1,:) = sum( reshape(GCOORD_SPH(3,EL2NOD(:,curved_els))',10,[]).*reshape(dN10(1,:,:),10,[]), 1);
        dTHdS(3,2,:) = sum( reshape(GCOORD_SPH(3,EL2NOD(:,curved_els))',10,[]).*reshape(dN10(2,:,:),10,[]), 1);
        dTHdS(3,3,:) = sum( reshape(GCOORD_SPH(3,EL2NOD(:,curved_els))',10,[]).*reshape(dN10(3,:,:),10,[]), 1);

        % determinants of the dTHdS matrices
        det_dTHdS = dTHdS(1,1,:).*dTHdS(2,2,:).*dTHdS(3,3,:) ...
                  + dTHdS(1,2,:).*dTHdS(2,3,:).*dTHdS(3,1,:) ...
                  + dTHdS(1,3,:).*dTHdS(2,1,:).*dTHdS(3,2,:) ...
                  - dTHdS(1,3,:).*dTHdS(2,2,:).*dTHdS(3,1,:) ...
                  - dTHdS(1,1,:).*dTHdS(2,3,:).*dTHdS(3,2,:) ...
                  - dTHdS(1,2,:).*dTHdS(2,1,:).*dTHdS(3,3,:);

        % 1/det(dXdSs)
        invdet = 1./det_dTHdS;

        inv_dTHdS(1,1,:) = invdet .* ( dTHdS(2,2,:).*dTHdS(3,3,:) - dTHdS(2,3,:).*dTHdS(3,2,:) );
        inv_dTHdS(2,1,:) = invdet .* ( dTHdS(2,3,:).*dTHdS(3,1,:) - dTHdS(2,1,:).*dTHdS(3,3,:) );
        inv_dTHdS(3,1,:) = invdet .* ( dTHdS(2,1,:).*dTHdS(3,2,:) - dTHdS(2,2,:).*dTHdS(3,1,:) );

        inv_dTHdS(1,2,:) = invdet .* ( dTHdS(1,3,:).*dTHdS(3,2,:) - dTHdS(1,2,:).*dTHdS(3,3,:) );
        inv_dTHdS(2,2,:) = invdet .* ( dTHdS(1,1,:).*dTHdS(3,3,:) - dTHdS(1,3,:).*dTHdS(3,1,:) );
        inv_dTHdS(3,2,:) = invdet .* ( dTHdS(1,2,:).*dTHdS(3,1,:) - dTHdS(1,1,:).*dTHdS(3,2,:) );

        inv_dTHdS(1,3,:) = invdet .* ( dTHdS(1,2,:).*dTHdS(2,3,:) - dTHdS(1,3,:).*dTHdS(2,2,:) );
        inv_dTHdS(2,3,:) = invdet .* ( dTHdS(1,3,:).*dTHdS(2,1,:) - dTHdS(1,1,:).*dTHdS(2,3,:) );
        inv_dTHdS(3,3,:) = invdet .* ( dTHdS(1,1,:).*dTHdS(2,2,:) - dTHdS(1,2,:).*dTHdS(2,1,:) );

        delS(1,:) =  inv_dTHdS(1,1,:) .* reshape(delTH(1,:),1,1,[]) ...
                   + inv_dTHdS(1,2,:) .* reshape(delTH(2,:),1,1,[]) ...
                   + inv_dTHdS(1,3,:) .* reshape(delTH(3,:),1,1,[]);
        delS(2,:) =  inv_dTHdS(2,1,:) .* reshape(delTH(1,:),1,1,[]) ...
                   + inv_dTHdS(2,2,:) .* reshape(delTH(2,:),1,1,[]) ...
                   + inv_dTHdS(2,3,:) .* reshape(delTH(3,:),1,1,[]);
        delS(3,:) =  inv_dTHdS(3,1,:) .* reshape(delTH(1,:),1,1,[]) ...
                   + inv_dTHdS(3,2,:) .* reshape(delTH(2,:),1,1,[]) ...
                   + inv_dTHdS(3,3,:) .* reshape(delTH(3,:),1,1,[]);
    end

    lc = lc + delS; % update guess
end

end % END OF FUNCTION local_coords_curved_3d_sph

% ##############################################################################################################################################################################################
%                               SUBFUNCTIONS
% ##############################################################################################################################################################################################

function [N10,dN10] = shapefun_10nodeVEC(r,s,t)
% Find shape functions and their derivatives at given points on the
% master element for a 10 node tetrahedron

% Node notation taken from Hughes' book, p.171;
% see also Zienkiewicz (4th) Volume 1, p.137
% local coords of the 10 nodes:
% Node  r  s  t
%  1    1  0  0
%  2    0  1  0
%  3    0  0  1
%  4    0  0  0
%  5   .5 .5  0
%  6    0 .5 .5
%  7    0  0 .5
%  8   .5  0  0
%  9   .5  0 .5
% 10    0 .5  0

npts = length(r);

N10  = zeros(10,npts);
dN10 = zeros(3,10,npts);
u    = 1-r-s-t;

% 4 vertex nodes
N10( 1,:) = r.*(2.*r-1);
N10( 2,:) = s.*(2.*s-1);
N10( 3,:) = t.*(2.*t-1);
N10( 4,:) = u.*(2.*u-1);
% 6 edge nodes
N10( 5,:) = 4.*r.*s;
N10( 6,:) = 4.*s.*t;
N10( 7,:) = 4.*t.*u;
N10( 8,:) = 4.*r.*u;
N10( 9,:) = 4.*r.*t;
N10(10,:) = 4.*s.*u;

% derivatives (3 for each node)
dN10(1,1,:) = 4*r-1; % dN1/dr
dN10(2,1,:) = 0;     % dN1/ds
dN10(3,1,:) = 0;     % dN1/dt

dN10(1,2,:) = 0;     % dN2/dr
dN10(2,2,:) = 4*s-1; % dN2/ds
dN10(3,2,:) = 0;     % dN2/dt

dN10(1,3,:) = 0;     % dN3/dr
dN10(2,3,:) = 0;     % dN3/ds
dN10(3,3,:) = 4*t-1; % dN3/dt

dN10(1,4,:) = -3+4*(1-u); % dN4/dr
dN10(2,4,:) = -3+4*(1-u); % dN4/ds
dN10(3,4,:) = -3+4*(1-u); % dN4/dt

dN10(1,5,:) = 4*s;   % dN5/dr
dN10(2,5,:) = 4*r;   % dN5/ds
dN10(3,5,:) = 0;     % dN5/dt

dN10(1,6,:) = 0;     % dN6/dr
dN10(2,6,:) = 4*t;   % dN6/ds
dN10(3,6,:) = 4*s;   % dN6/dt

dN10(1,7,:) = -4*t;    % dN7/dr
dN10(2,7,:) = -4*t;    % dN7/ds
dN10(3,7,:) = 4*(u-t); % dN7/dt

dN10(1,8,:) = 4*(u-r); % dN8/dr
dN10(2,8,:) = -4*r;    % dN8/ds
dN10(3,8,:) = -4*r;    % dN8/dt

dN10(1,9,:) = 4*t;   % dN9/dr
dN10(2,9,:) = 0;     % dN9/ds
dN10(3,9,:) = 4*r;   % dN9/dt

dN10(1,10,:) = -4*s;    % dN10/dr
dN10(2,10,:) = 4*(u-s); % dN10/ds
dN10(3,10,:) = -4*s;    % dN10/dt

end % END OF FUNCTION shapefun_10nodeVEC