function RR = make_rot_matrix_freeslip(MESH,PlateInfo,useNataf)
%
% Purpose:
% Create the transformation(=all local rotations) matrix RR that transforms
% (x,y,z) coordinate systems for points on Core-Mantle Boundary (CMB) to 
% local (r,N,E) coordinate system.  This is needed to create a 
% shear-stress free BC at the CMB, where only vR=0 will be a prescribed BC.
%
% Input: 
%   MESH :: structure containing all FE MESH information
%           used for its GCOORD locations of node points and their DBnods values
% Output:
%   RR :: transformation matrix that rotates all global dofs to their preferred
%         local coordinate system
% JPM April 2011

GCOORD   = MESH.GCOORD;
PointID  = MESH.PointID;
nnod     = MESH.nnod;
ndof     = 3;
nUdof    = nnod * ndof;
nods_CMB = find(PointID==301);

if useNataf
    % Free slip nodes IF asth BC is assumed
    nods_TopFreeSlip = PlateInfo(1, PlateInfo(2,:)>=3 );
else
    nods_TopFreeSlip = []; % DONT ADD ANY TOP FREE SLIP nodes
end

nods_CMB = [nods_CMB; nods_TopFreeSlip];

nCMBnod = length(nods_CMB);              % # of CMB nodes
scale   = sum( GCOORD(:,nods_CMB).^2,1); % dot product -- unrolled loop for l^2
% unit vectors at locations of CMB nodes
% scale by rCMB to make unit vector rhat
rhatCMB = (GCOORD(:,nods_CMB) ./ repmat( sqrt(scale) ,[3 1])); 

% Turn all rhatCMB unit vectors into latlon
[xyzlatlon] = vangV(rhatCMB); % *SUBFUNCTION*

% now make vectors xyzN that are 90deg N of each point (wrapping over the pole OK for this logic)
[xyzN] = vsetV(xyzlatlon(1,:)+90.0,xyzlatlon(2,:)); % *SUBFUNCTION*

% now make vectors xyzE pointing east at each point (compute at lat=equator for easier calculation)
[xyzE] = vsetV(zeros(1,size(xyzlatlon,2)),xyzlatlon(2,:)+90.0); % *SUBFUNCTION* 

% rNE coordinate unit vectors for each point are rhatCMB, xyzN, xyzE
% these define the rotation matrix components for each point on CMB
RR_CMB = zeros(3,3,nCMBnod);
RR_CMB(1,1,:) = rhatCMB(1,:); RR_CMB(1,2,:) = rhatCMB(2,:);RR_CMB(1,3,:) = rhatCMB(3,:);
RR_CMB(2,1,:) = xyzN(1,:);    RR_CMB(2,2,:) = xyzN(2,:);   RR_CMB(2,3,:) = xyzN(3,:);
RR_CMB(3,1,:) = xyzE(1,:);    RR_CMB(3,2,:) = xyzE(2,:);   RR_CMB(3,3,:) = xyzE(3,:);

% now need to put these RR_CMB matrices into the proper global dof
% locations in a global transformation matrix RR that transforms all
% local coordinate directions (either (r,N,E) at CMB or (x,y,z) elsewhere)
% into global ones ((x,y,z) everywhere). For all nonCMB dofs, this will
% just be an identity matrix, for CMB nodes it will be the corresponding
% entries in RR_CMB.

% pre-allocate number of non-zero slots in RR (RR is diagonal with three
% dofs for each node, EXCEPT for 3x3 submatrices for each of the nodes in CMB)
nnz_RR = nUdof + 6*nCMBnod; % 6 extra dof needed to fill in 3x3 submatrices for CMB pts
RRi = zeros(nnz_RR,1); % initialize nonzero i,j, value vectors 
RRj = zeros(nnz_RR,1); 
RRv = ones(nnz_RR,1);  % NOTE  initialize to possible Identity matrix value 

%on-diagonal entries of RR: This will use up the first 3*npt entries
% identity matrix part of entries of T (RRv already set to 1)
RRi(1:nUdof) = [1:nUdof]';
RRj(1:nUdof) = [1:nUdof]';

% ---
% RR_CMB = ones(size(RR_CMB)); % use 1) this + 2) spy command to check sparsity pattern
% ---

% Now overwrite RR_CMB's onto the diagonal entries for the CMB points
% CMBdofs = [3*(nods_CMB-1)+1;3*(nods_CMB-1)+2;3*(nods_CMB-1)+3];
RRv(3*(nods_CMB-1)+1) = RR_CMB(1,1,:); %x'-entry diagonal
RRv(3*(nods_CMB-1)+2) = RR_CMB(2,2,:); %y'-entry diagonal
RRv(3*(nods_CMB-1)+3) = RR_CMB(3,3,:); %z'-entry diagonal



% Now add the 6 off-diagonal slots for CMB points
nbeg = nUdof+1;
nend = nbeg+nCMBnod-1; 
RRi(nbeg:nend) = [3*(nods_CMB-1)+1];   % i for 1,2 dofs entry
RRj(nbeg:nend) = [3*(nods_CMB-1)+2];   % j for 1,2 dofs entry  
RRv(nbeg:nend) = RR_CMB(1,2,:);       % value 

nbeg = nend+1;
nend = nbeg+nCMBnod-1; 
RRi(nbeg:nend) = [3*(nods_CMB-1)+1];   % i for 1,3 dofs entry
RRj(nbeg:nend) = [3*(nods_CMB-1)+3];   % j for 1,3 dofs entry  
RRv(nbeg:nend) = RR_CMB(1,3,:);       % value 

nbeg = nend+1;
nend = nbeg+nCMBnod-1; 
RRi(nbeg:nend) = [3*(nods_CMB-1)+2];   % i for 2,3 dofs entry
RRj(nbeg:nend) = [3*(nods_CMB-1)+3];   % j for 2,3 dofs entry  
RRv(nbeg:nend) = RR_CMB(2,3,:);       % value 

nbeg = nend+1;
nend = nbeg+nCMBnod-1; 
RRi(nbeg:nend) = [3*(nods_CMB-1)+2];   % i for 2,1 dofs entry
RRj(nbeg:nend) = [3*(nods_CMB-1)+1];   % j for 2,1 dofs entry  
RRv(nbeg:nend) = RR_CMB(2,1,:);       % value 

nbeg = nend+1;
nend = nbeg+nCMBnod-1; 
RRi(nbeg:nend) = [3*(nods_CMB-1)+3];   % i for 3,1 dofs entry
RRj(nbeg:nend) = [3*(nods_CMB-1)+1];   % j for 3,1 dofs entry  
RRv(nbeg:nend) = RR_CMB(3,1,:);       % value 

nbeg = nend+1;
nend = nbeg+nCMBnod-1; 
RRi(nbeg:nend) = [3*(nods_CMB-1)+3];   % i for 3,2 dofs entry
RRj(nbeg:nend) = [3*(nods_CMB-1)+2];   % j for 3,2 dofs entry  
RRv(nbeg:nend) = RR_CMB(3,2,:);       % value 

RR = sparse(RRi(:),RRj(:),RRv(:)); % MATLAB sparse assembles from i,j,v triplets

end % END OF FUNCTION make_rot_matrix

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function [alatlong] = vangV(a)
%c        convert from direction cosines into  lat,long   26apr99
%c        (vector 'a' does not have to be of unit length.)
  rad2deg =  1./0.0174532925199;
  aa    = zeros(1,size(a,2));   % initialize variables
  alat  = zeros(1,size(a,2));
  along = zeros(1,size(a,2));
  
  aa(:)    = sqrt( a(1,:).*a(1,:) + a(2,:).*a(2,:) );
  alat(:)  = 90.0 .* sign( a(3,:) );
  along(:) = 0.0;
  
  alat(aa > 0.0)  = atan ( a(3,aa > 0.0) ./ aa(aa > 0.0)  ) .* rad2deg;
  along(aa > 0.0) = atan2( a(2,aa > 0.0), a(1,aa > 0.0) ) .* rad2deg;

%c...(uncomment next lines to go from 0-360 degrees instead of -180 to +180)
  along(along < 0.0 )     = along(along < 0.0 ) + 360.;
  along(along > 359.9999) = 0.;
  alatlong = [alat;along];
end  % end of function vangV

% ###################################################################

function [a] = vsetV(alat,along)
%c convert  lat,long in degrees into direction cosines  26apr99
  a = zeros(3,max(size(alat)));
  deg2rad =  0.0174532925199;
  a(1,:) = cos(deg2rad*alat(:)) .* cos(deg2rad*along(:));
  a(2,:) = cos(deg2rad*alat(:)) .* sin(deg2rad*along(:));
  a(3,:) = sin(deg2rad*alat(:));
end % end of function vset