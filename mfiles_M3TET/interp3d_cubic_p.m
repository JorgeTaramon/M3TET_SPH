function [V_i,OPTS] = interp3d_cubic_p(GCOORD,EL2NOD,COMM,els,lc,V,OPTS)
% Usage: [V_i,OPTS] = interp3d_cubic_p(GCOORD,EL2NOD,COMM,els,lc,V,OPTS)
%
% Purpose: The function performs a quasi-cubic interpolation on an unstructured
%          triangular mesh. Based on code developed by Chao Shi & Jason
%          Phipps Morgan, 2008.
% Input:
%   GCOORD : [matrix] : coordinates of all nodes in mesh (3 x nnodel)
%   EL2NOD : [matrix] : finite element connectivity matrix (nnodel x nel)
%   COMM   : [structure] : inter-subdomain communication data
%   els    : [vector] : element in which each point is located (npt x 1)
%   lc     : [matrix] : local coordinates in each element (3 x npt)
%   V      : [matrix] : variable field(s) to be interpolated (nnod x nvar)
%   OPTS   : [structure] : options for this functions
%
% Output
%   V_i    : [matrix] : interpolated values (npt x nvar)
%   OPTS   : [structure] : options for this functions
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Written by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)

% JH April 2011
% JH March 2013
% JH August 2013
% JH April 2015
%

if nargin==6; OPTS = []; end

% DEFAULT SETTINGS (will be overwritten if options (OPTS) are provided
if ~isfield(OPTS,'verbose')
    OPTS.verbose = 0; % display profiling data
end
if ~isfield(OPTS,'monotonic')
    OPTS.monotonic = 1; % cut-off for interpolated values
end
if ~isfield(OPTS,'method_wght')
    % 'ones' :: uniform weighting
    % 'dist' :: weighted by distance of node to element center
    OPTS.method_wght = 'dist'; 
end
if ~isfield(OPTS,'nelblk')
    OPTS.nelblk = 25000;  % number of elements in block
end


nnod = size(GCOORD,2);    % number of nodes in mesh
if size(V,2)==nnod; V=V'; end
nvar = size(V,2);
els  = els(:);

% nelblk  total time (L190/ 321 / 360)
%  1000 --> 26.7 sec (1.0 / 0.25/ 5.3)
%  5000 --> 10.5 sec (1.0 / 0.3 / 1.0) % best nelblk for invM*
% 10000 --> 10.5 sec (1.2 / 0.4 / 0.5)
% 20000 --> 11.5 sec (2.2 / 0.5 / 0.2)
% 30000 --> 10.6 sec (2.3 / 0.5 / 0.2)
% 50000 --> 10.0 sec (2.5 / 0.4 / 0.1) % best nelblk for dBdxyz_nod(:,*) + tmp

                                                                            tic
% Calculates the following spacial derivatives:
%  dV1  dV1  dV1
% ---- ---- ----
%  dx   dy   dz

%NOTICE: the reason for adding sqrt(2) is in the notes
%==========================================================================
% Shape functions 5,6,9,11,12,15 are different b/c they are on tilted sides
% We add sqrt(2) to proper fit the slopeDev information, see notes for the
% math
%==========================================================================
%
%                                      -- Chao 12/03/2008
%
r = lc(1,:);
s = lc(2,:);
t = lc(3,:);

% shape functions for cubic interpolation
npt = length(els);
if npt
    N  = zeros(16,npt);
    u  = 1-r-s-t;
    % linear
    N(1,:)  = r;
    N(2,:)  = s;
    N(3,:)  = t;
    N(4,:)  = u;
    % quadratic
    N(5,:)  = sqrt(2)*r.*s;
    N(6,:)  = sqrt(2)*s.*t;
    N(7,:)  = t.*u;
    N(8,:)  = r.*u;
    N(9,:)  = sqrt(2)*r.*t;
    N(10,:) = s.*u;
    % cubic
    N(11,:) = sqrt(2)*s.*r.*(s-r);
    N(12,:) = sqrt(2)*s.*t.*(t-s);
    N(13,:) = t.*u.*(u-t);
    N(14,:) = r.*u.*(u-r);
    N(15,:) = sqrt(2)*r.*t.*(t-r);
    N(16,:) = s.*u.*(u-s);
%     % bubble (not used here)
%     N(17,:) = 27 * s.*t.*u;
%     N(18,:) = 27 * r.*t.*u;
%     N(19,:) = 27 * r.*s.*u;
%     N(20,:) = 27 * r.*s.*t;
else
    N = [];
end

% Nodal coordinates; four x values for a tetrahedron in each column, etc
x = reshape(GCOORD(1,EL2NOD(1:4,els)),4,npt);
y = reshape(GCOORD(2,EL2NOD(1:4,els)),4,npt);
z = reshape(GCOORD(3,EL2NOD(1:4,els)),4,npt);

return_wght = 1;
V_i = zeros(npt,nvar);
for ivar=1:nvar % loop over variable fields
    B = V(:,ivar);
    
    % Calculate spatial derivatives of variable
    [dBdxyz,dBdxyz_el,tmp] = calc_derivatives...
        (GCOORD,EL2NOD(1:4,:),B,OPTS.nelblk,OPTS.method_wght,return_wght); % *SUBFUNCTION*
    if return_wght; wght = tmp; end
    
    % PARALLEL SECTION: Communication required to calculate spatial
    % derivatives of "B" and the weights (only once) at SD boundaries
    if COMM.nsd>1
        if return_wght; wght_MY = wght; end % Must be copied because it will be updated within the
        dBdxyz_MY = dBdxyz;                 % ...loop below but original values must be sent to NBs
        mynod_SDB = COMM.mynod_SDB{1};
        for iNB=1:COMM.nNB
            labBarrier
            NB = COMM.NB(iNB);
            if NB>0
                SDBnod = mynod_SDB{iNB};
                if return_wght
                    wght_NB      = COMM.sendnrecv_1NB(NB,wght_MY(SDBnod),201);
                    wght(SDBnod) = wght(SDBnod) + wght_NB; clear wght_NB
                end
                dBdxyz_NB        = COMM.sendnrecv_1NB(NB,dBdxyz_MY(SDBnod,:),202);
                dBdxyz(SDBnod,:) = dBdxyz(SDBnod,:) + dBdxyz_NB; clear dBdxyz_NB
            end
        end
        clear dBdxyz_MY dBdxyz_NB
    end
    return_wght = 0;

    % NOTE: Even if "els" is empty (no points to be interpolated in this SD) 
    %       the parallel code must run to this point because of the
    %       communication above!! However, the calculations below can 
    %       be skipped.
    if isempty(els)
        continue
    end
    
    dBdxyz = dBdxyz ./ repmat(wght,1,3);
    dBdx   = reshape(dBdxyz(EL2NOD(1:4,els),1),4,npt) - repmat(dBdxyz_el(els,1)',4,1);
    dBdy   = reshape(dBdxyz(EL2NOD(1:4,els),2),4,npt) - repmat(dBdxyz_el(els,2)',4,1);
    dBdz   = reshape(dBdxyz(EL2NOD(1:4,els),3),4,npt) - repmat(dBdxyz_el(els,3)',4,1);

    % Pre-allocate memory
    slopeDev_r = zeros(4,npt);
    slopeDev_s = zeros(4,npt);
    slopeDev_t = zeros(4,npt);

    % Map the slope DEVIATION from (x,y,z) space to (r,s,t) space
    % Jacobian matrix: J = [1 0 0 -1;0 1 0 -1;0 0 1 -1] * x_y_z matrix
    slopeDev_r(1,:) = (x(1,:)-x(4,:)).*dBdx(1,:)+(y(1,:)-y(4,:)).*dBdy(1,:)+(z(1,:)-z(4,:)).*dBdz(1,:);
    slopeDev_r(2,:) = (x(1,:)-x(4,:)).*dBdx(2,:)+(y(1,:)-y(4,:)).*dBdy(2,:)+(z(1,:)-z(4,:)).*dBdz(2,:);
    slopeDev_r(3,:) = (x(1,:)-x(4,:)).*dBdx(3,:)+(y(1,:)-y(4,:)).*dBdy(3,:)+(z(1,:)-z(4,:)).*dBdz(3,:);
    slopeDev_r(4,:) = (x(1,:)-x(4,:)).*dBdx(4,:)+(y(1,:)-y(4,:)).*dBdy(4,:)+(z(1,:)-z(4,:)).*dBdz(4,:);

    slopeDev_s(1,:) = (x(2,:)-x(4,:)).*dBdx(1,:)+(y(2,:)-y(4,:)).*dBdy(1,:)+(z(2,:)-z(4,:)).*dBdz(1,:);
    slopeDev_s(2,:) = (x(2,:)-x(4,:)).*dBdx(2,:)+(y(2,:)-y(4,:)).*dBdy(2,:)+(z(2,:)-z(4,:)).*dBdz(2,:);
    slopeDev_s(3,:) = (x(2,:)-x(4,:)).*dBdx(3,:)+(y(2,:)-y(4,:)).*dBdy(3,:)+(z(2,:)-z(4,:)).*dBdz(3,:);
    slopeDev_s(4,:) = (x(2,:)-x(4,:)).*dBdx(4,:)+(y(2,:)-y(4,:)).*dBdy(4,:)+(z(2,:)-z(4,:)).*dBdz(4,:);

    slopeDev_t(1,:) = (x(3,:)-x(4,:)).*dBdx(1,:)+(y(3,:)-y(4,:)).*dBdy(1,:)+(z(3,:)-z(4,:)).*dBdz(1,:);
    slopeDev_t(2,:) = (x(3,:)-x(4,:)).*dBdx(2,:)+(y(3,:)-y(4,:)).*dBdy(2,:)+(z(3,:)-z(4,:)).*dBdz(2,:);
    slopeDev_t(3,:) = (x(3,:)-x(4,:)).*dBdx(3,:)+(y(3,:)-y(4,:)).*dBdy(3,:)+(z(3,:)-z(4,:)).*dBdz(3,:);
    slopeDev_t(4,:) = (x(3,:)-x(4,:)).*dBdx(4,:)+(y(3,:)-y(4,:)).*dBdy(4,:)+(z(3,:)-z(4,:)).*dBdz(4,:);

    % For edges not parallel to any of the axis (r,s,t), we need to find
    % the slope along the tilted edges
    dTdq1node2 = sqrt(1/2)*slopeDev_r(2,:) - sqrt(1/2)*slopeDev_s(2,:);
    dTdq1node1 = sqrt(1/2)*slopeDev_r(1,:) - sqrt(1/2)*slopeDev_s(1,:);

    dTdq2node3 = sqrt(1/2)*slopeDev_s(3,:) - sqrt(1/2)*slopeDev_t(3,:);
    dTdq2node2 = sqrt(1/2)*slopeDev_s(2,:) - sqrt(1/2)*slopeDev_t(2,:);

    dTdq5node3 = sqrt(1/2)*slopeDev_r(3,:) - sqrt(1/2)*slopeDev_t(3,:);
    dTdq5node1 = sqrt(1/2)*slopeDev_r(1,:) - sqrt(1/2)*slopeDev_t(1,:);

    % Calculate the slope information for each edge
    slopeDev_info       = zeros(12,npt);
    slopeDev_info(1,:)  = (dTdq1node2-dTdq1node1)/2;              % for N5
    slopeDev_info(2,:)  = (dTdq2node3-dTdq2node2)/2;              % for N6
    slopeDev_info(3,:)  = (slopeDev_t(4,:)-slopeDev_t(3,:))/2;    % for N7
    slopeDev_info(4,:)  = (slopeDev_r(4,:)-slopeDev_r(1,:))/2;    % for N8
    slopeDev_info(5,:)  = (dTdq5node3-dTdq5node1)/2;              % for N9
    slopeDev_info(6,:)  = (slopeDev_s(4,:)-slopeDev_s(2,:))/2;    % for N10
    slopeDev_info(7,:)  = (dTdq1node2+dTdq1node1)/2;              % for N11
    slopeDev_info(8,:)  = (dTdq2node3+dTdq2node2)/2;              % for N12
    slopeDev_info(9,:)  = (slopeDev_t(4,:)+slopeDev_t(3,:))/2;    % for N13
    slopeDev_info(10,:) = (slopeDev_r(4,:)+slopeDev_r(1,:))/2;    % for N14
    slopeDev_info(11,:) = (dTdq5node3+dTdq5node1)/2;              % for N15
    slopeDev_info(12,:) = (slopeDev_s(4,:)+slopeDev_s(2,:))/2;    % for N16

    B_i = sum( N(1:16,:).*[B(EL2NOD(1:4,els)) ; slopeDev_info(1:12,:)] )';

    if OPTS.monotonic
        B_elnods = B(EL2NOD(1:4,:));
        B_max_el = max(B_elnods,[],1)';
        B_min_el = min(B_elnods,[],1)';
        B_max    = B_max_el(els);
        B_min    = B_min_el(els);
        B_i      = max(B_i,B_min);
        B_i      = min(B_i,B_max);
    end
    
    V_i(:,ivar) = B_i;
end

if OPTS.verbose
    fprintf(' Cubic interpolation (monotonic=%1i,nelblk=%1i,wght=%s):\n',...
            OPTS.monotonic,OPTS.nelblk,OPTS.method_wght);
    fprintf('      %1i variable(s) at %1i points in %5.2e sec\n',...
                 nvar,npt,toc);
end

end % END OF FUNCTION interp3d_cubic_p

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function [dBdxyz_nod,dBdxyz_el,wght_nod] = calc_derivatives...
             (GCOORD,EL2NOD,B,nelblk,method_wght,return_wght)

% The following lines are based on the performance-improved assembly described 
% in the paper "MILAMIN: Matlab-based FEM solver for large problems" by
% Dabrowski, Krotkiewski and Schmid; G3, VOL. 9, doi:10.1029/2007GC001719, 2008
nel    = size(EL2NOD,2);
nnod   = size(GCOORD,2);   % number of nodes in mesh
nelblk = min(nel,nelblk);  % in case blocksize>number of elements
nblk   = ceil(nel/nelblk); % number of blocks

invM1  = zeros(nelblk,4);  % invM1 stores the 1st row of inv(M) of all elements in block
invMx  = zeros(nelblk,4);  % invMx stores the 2nd row...
invMy  = zeros(nelblk,4);  % invMy stores the 3rd row...
invMz  = zeros(nelblk,4);  % invMz stores the 4th row...

dBdxyz_el  = zeros(nel,3);
dBdxyz_nod = zeros(nnod,3);
if return_wght
    wght_nod = zeros(nnod,1);
else
    wght_nod = [];
end

il = 1;       % 1st element in 1st block
iu = nelblk;  % last element in 1st block
for ib=1:nblk % loop over the element blocks
    
    EL2NOD_BLK = EL2NOD(:,il:iu)';
    
    M1 = ones(nelblk,4);       % First column of "M" for all elements in block
    Mx = GCOORD(1,EL2NOD_BLK); % 2nd column of "M"... (= x-coordinates)
    Mx = reshape(Mx,nelblk,4); % 1st row = 1st element in block; 2nd row...
    My = GCOORD(2,EL2NOD_BLK); % 3rd column of "M"... (= y-coordinates)
    My = reshape(My,nelblk,4); % 1st row = 1st element in block; 2nd row...
    Mz = GCOORD(3,EL2NOD_BLK); % 4th column of "M"... (= z-coordinates)
    Mz = reshape(Mz,nelblk,4); % 1st row = 1st element in block; 2nd row...
    
    switch method_wght
        case 'ones'
            wght = ones(nelblk,4);
        case 'dist'
            Cx   = repmat(sum(Mx,2)./4,1,4); % x coord of element centers
            Cy   = repmat(sum(My,2)./4,1,4); % y coord of element centers
            Cz   = repmat(sum(Mz,2)./4,1,4); % z coord of element centers
            wght = sqrt( (Mx-Cx).^2 + (My-Cy).^2 + (Mz-Cz).^2);
            wght = 1./wght;
    end
    
    % Determinate of "M" for all elements in block, dropped M1(all ones)
    detM = Mx(:,2).*My(:,3).*Mz(:,4) - Mx(:,2).*Mz(:,3).*My(:,4)...
         - Mx(:,3).*My(:,2).*Mz(:,4) + Mx(:,3).*Mz(:,2).*My(:,4)...
         + Mx(:,4).*My(:,2).*Mz(:,3) - Mx(:,4).*Mz(:,2).*My(:,3)...
         - Mx(:,1).*My(:,3).*Mz(:,4) + Mx(:,1).*Mz(:,3).*My(:,4)...
         + Mx(:,3).*My(:,1).*Mz(:,4) - Mx(:,3).*Mz(:,1).*My(:,4)...
         - Mx(:,4).*My(:,1).*Mz(:,3) + Mx(:,4).*Mz(:,1).*My(:,3)...
         + Mx(:,1).*My(:,2).*Mz(:,4) - Mx(:,1).*Mz(:,2).*My(:,4)...
         - Mx(:,2).*My(:,1).*Mz(:,4) + Mx(:,2).*Mz(:,1).*My(:,4)...
         + Mx(:,4).*My(:,1).*Mz(:,2) - Mx(:,4).*Mz(:,1).*My(:,2)...
         - Mx(:,1).*My(:,2).*Mz(:,3) + Mx(:,1).*Mz(:,2).*My(:,3)...
         + Mx(:,2).*My(:,1).*Mz(:,3) - Mx(:,2).*Mz(:,1).*My(:,3)...
         - Mx(:,3).*My(:,1).*Mz(:,2) + Mx(:,3).*Mz(:,1).*My(:,2);
     
	% 1/det(M) for all elements in block
    invdetM = 1./detM;
     
    % Inversion of M matrices for all elements in block
    % e.g.: invM1 is 1st row of inv(M) of all elements in block
    invM1(:,1) =   invdetM .* det33b(Mx(:,2),Mx(:,3),Mx(:,4),My(:,2),My(:,3),My(:,4),Mz(:,2),Mz(:,3),Mz(:,4));
    invM1(:,2) = - invdetM .* det33b(Mx(:,1),Mx(:,3),Mx(:,4),My(:,1),My(:,3),My(:,4),Mz(:,1),Mz(:,3),Mz(:,4));
    invM1(:,3) =   invdetM .* det33b(Mx(:,1),Mx(:,2),Mx(:,4),My(:,1),My(:,2),My(:,4),Mz(:,1),Mz(:,2),Mz(:,4));
    invM1(:,4) = - invdetM .* det33b(Mx(:,1),Mx(:,2),Mx(:,3),My(:,1),My(:,2),My(:,3),Mz(:,1),Mz(:,2),Mz(:,3));

    invMx(:,1) = - invdetM .* det33b(M1(:,2),M1(:,3),M1(:,4),My(:,2),My(:,3),My(:,4),Mz(:,2),Mz(:,3),Mz(:,4));
    invMx(:,2) =   invdetM .* det33b(M1(:,1),M1(:,3),M1(:,4),My(:,1),My(:,3),My(:,4),Mz(:,1),Mz(:,3),Mz(:,4));
    invMx(:,3) = - invdetM .* det33b(M1(:,1),M1(:,2),M1(:,4),My(:,1),My(:,2),My(:,4),Mz(:,1),Mz(:,2),Mz(:,4));
    invMx(:,4) =   invdetM .* det33b(M1(:,1),M1(:,2),M1(:,3),My(:,1),My(:,2),My(:,3),Mz(:,1),Mz(:,2),Mz(:,3));

    invMy(:,1) =   invdetM .* det33b(M1(:,2),M1(:,3),M1(:,4),Mx(:,2),Mx(:,3),Mx(:,4),Mz(:,2),Mz(:,3),Mz(:,4));
    invMy(:,2) = - invdetM .* det33b(M1(:,1),M1(:,3),M1(:,4),Mx(:,1),Mx(:,3),Mx(:,4),Mz(:,1),Mz(:,3),Mz(:,4));
    invMy(:,3) =   invdetM .* det33b(M1(:,1),M1(:,2),M1(:,4),Mx(:,1),Mx(:,2),Mx(:,4),Mz(:,1),Mz(:,2),Mz(:,4));
    invMy(:,4) = - invdetM .* det33b(M1(:,1),M1(:,2),M1(:,3),Mx(:,1),Mx(:,2),Mx(:,3),Mz(:,1),Mz(:,2),Mz(:,3));

    invMz(:,1) = - invdetM .* det33b(M1(:,2),M1(:,3),M1(:,4),Mx(:,2),Mx(:,3),Mx(:,4),My(:,2),My(:,3),My(:,4));
    invMz(:,2) =   invdetM .* det33b(M1(:,1),M1(:,3),M1(:,4),Mx(:,1),Mx(:,3),Mx(:,4),My(:,1),My(:,3),My(:,4));
    invMz(:,3) = - invdetM .* det33b(M1(:,1),M1(:,2),M1(:,4),Mx(:,1),Mx(:,2),Mx(:,4),My(:,1),My(:,2),My(:,4));
    invMz(:,4) =   invdetM .* det33b(M1(:,1),M1(:,2),M1(:,3),Mx(:,1),Mx(:,2),Mx(:,3),My(:,1),My(:,2),My(:,3));

    
    % B_blk is the variable value at the 3 nodes of all elements in block
    % ==> has same dimensions as invM1, invMx, invMy, invMz: nelblk x 3
    % Thus, we can perform a point-wise multiplication .* and sum the rows
    
    % Variable value at 3 nodes of all elements in block;
    % reshape it so this line works for only 1 element
    B_blk    = reshape(B(EL2NOD_BLK),nelblk,[]);
    dBdx_blk = sum( invMx.*B_blk , 2 ); % 2nd coefficient == invM(2,:) .* nodal values
    dBdy_blk = sum( invMy.*B_blk , 2 ); % 3rd coefficient == invM(3,:) .* nodal values
    dBdz_blk = sum( invMz.*B_blk , 2 ); % 4th coefficient == invM(4,:) .* nodal values

    tmp1            = wght.*repmat(dBdx_blk,1,4);
    tmp2            = accumarray( EL2NOD_BLK(:), tmp1(:), [nnod 1] );
    dBdxyz_nod(:,1) = dBdxyz_nod(:,1) + tmp2;
    
    tmp1            = wght.*repmat(dBdy_blk,1,4);
    tmp2            = accumarray( EL2NOD_BLK(:), tmp1(:), [nnod 1] );
    dBdxyz_nod(:,2) = dBdxyz_nod(:,2) + tmp2;
    
    tmp1            = wght.*repmat(dBdz_blk,1,4);
    tmp2            = accumarray( EL2NOD_BLK(:), tmp1(:), [nnod 1] );
    dBdxyz_nod(:,3) = dBdxyz_nod(:,3) + tmp2;
    
    dBdxyz_el(il:iu,:) = [dBdx_blk dBdy_blk dBdz_blk];
    
    if return_wght
        tmp2     = accumarray( EL2NOD_BLK(:), wght(:), [nnod 1]);
        wght_nod = wght_nod + tmp2;
    end
    
    % Go to next block of elements (increase element indices)
    il = il+nelblk;
    if (ib==nblk-1) % last block has to be scaled to the number of elements remaining
        nelblk = nel-iu;
        invM1  = zeros(nelblk,4);
        invMx  = zeros(nelblk,4);
        invMy  = zeros(nelblk,4);
        invMz  = zeros(nelblk,4);
    end
    iu = iu + nelblk;
end

end % END OF SUBFUNCTION calc_derivatives

% #########################################################################

function det_3by3 = det33b(a11,a12,a13,a21,a22,a23,a31,a32,a33)
    det_3by3 = a11.*a22.*a33 - a11.*a23.*a32 - a21.*a12.*a33...
             + a21.*a13.*a32 + a31.*a12.*a23 - a31.*a13.*a22;
end % END OF SUBFUNCTION det33b