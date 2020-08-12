function Ic2f = intermesh_transfer_matrix(GCOORD,EL2NOD,GCOORD_c,EL2NOD_c,ndof)

% (1) Locate nodes of fine mesh in coarse mesh elements
% =====================================================
% Calculate coarse element within which each fine node is located
% This does not require a call of tsearchn since the fine mesh is nested in
% the coarse one.
[nod_f,loc] = unique(EL2NOD);
% "nod_f" is a vector from 1:nnod.
% "loc" is the linear index (see "help sub2ind") corresponding to the last
% row/column entry in nodes with number inod (see below). The above line is
% equivalent to the following loop over all nodes on the fine mesh.
% loc = zeros(1,nnod);
% for inod=1:nnod
%     loc(inod) = find(EL2NOD==inod,1,'last');
% end

% Get the element number from loc
nnodel = size(EL2NOD,1);
els    = ceil(loc/nnodel);

% Calculate coarse element in which the fine element in "els" is nested.
els_c = ceil(els./8);

% (2) get local coordinate of each fine node in its coarse element
% ================================================================
lc = local_coords_3d(GCOORD_c,EL2NOD_c,els_c,GCOORD);
% Check if local coordinates are within the element
tol    = 1e-5;
if min(min(lc)) < -tol || max(max(lc)) > 1+tol || max(sum(lc)) > 1+tol
    error(' Calculation of interpolation matrix failed.');
end

% For a 10-node tetrahedra mesh, lc must only contain 5 different
% values: 0, 0.25, 0.5, 0.75, and 1
% To avoid unnecessary round-off errors in Ic2f, lc is adjusted to
% contain these 5 values only:
lc( lc<0 ) = 0;
lc( lc>1 ) = 1;
lc         = round(1000*lc)./1000;
% NOTE: If you don't cut of the round-off here, the restricted residual has
%       serial vs parallel errors of magnitude 1e-2 !!!!!!!

% Assembly for 1 d.o.f. per node (e.g. pressure in HT3-codes)
nnod_f = size(lc,2);
nnod_c = size(GCOORD_c,2);
ii     = repmat(nod_f,1,nnodel); % i-index in interpolation matrix
jj     = EL2NOD_c(:,els_c)';     % j-index in interpolation matrix
N      = sf_dsf_tet(lc,nnodel,'matrix')';
if ndof==3 % 3 dof per node
    ii     = nod2dof(ii(:),3);
    jj     = nod2dof(jj(:),3);
    N      = [N(:)'; N(:)'; N(:)']; N = N(:);
    ndof_f = 3*nnod_f;
    ndof_c = 3*nnod_c;
elseif ndof==1 % 1 dof per node
    ndof_f = nnod_f;
    ndof_c = nnod_c;
end

% Assemble triplet to form sparse matrix
Ic2f = sparse3(ii(:),jj(:),N(:),ndof_f,ndof_c);

% % Use this block for debugging purposes
% Ic2f_U_1 = zeros(nnod_f,nnod_c);
% Ic2f_U_3 = zeros(3*nnod_f,3*nnod_c);
% for inod_f=1:nnod_f
%     el_c                    = els_c(inod_f);
%     nods_c                  = EL2NOD_c(el_c,:);
%     Ic2f_U_1(inod_f,nods_c) = N(inod_f,:);
%     for idof=1:3
%         Ic2f_U_3(3*(inod_f-1)+idof,3*(nods_c-1)+idof) = N(inod_f,:);
%     end
% end
% varIc = GCOORD_c(1,:).*GCOORD_c(2,:);
% varIf = Ic2f_U_1 * varIc';
% max( abs((varIf-(GCOORD(1,:).*GCOORD(2,:))' )./(GCOORD(1,:).*GCOORD(2,:))') )

end % END OF SUBFUNCTION intermesh_transfer_matrix