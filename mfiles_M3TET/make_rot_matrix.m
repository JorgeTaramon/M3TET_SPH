function [RR,xyzR,xyzE,xyzN] = make_rot_matrix(GCOORD) %,PlateInfo,useNataf)
%
% Purpose:
% Create the transformation(=all local rotations) matrix RR that transforms
% (x,y,z) coordinate systems to local (r,E,N) coordinate system. This is 
% needed to create a shear-stress free BC at the CMB, where only vR=0 will 
% be a prescribed BC.
%
% Input: 
%   GCOORD :: xyz-coordinates of all nodes
%
% Output:
%   RR     :: transformation matrix that rotates all global dofs to their 
%             preferred local coordinate system
%
% JPM April 2011
% JH Mar 2015 : Now assembles rotation matrix for ALL dofs (not only dofs 
%               with free slip boundary condition)

nnod   = size(GCOORD,2);
ndof   = 3*nnod;
nods   = 1:nnod;

scale = sum(GCOORD(:,nods).^2,1);
% unit vectors at locations of nodes scale by radius to make unit vector rhat
rhat  = (GCOORD(:,nods) ./ repmat( sqrt(scale) ,[3 1])); 
xyzR  = rhat;

% Turn all rhat unit vectors into latlon
[xyzlatlon] = vangV(rhat); % *SUBFUNCTION*

% now make vectors xyzN that are 90deg N of each point
% (wrapping over the pole OK for this logic)
[xyzN] = vsetV(xyzlatlon(1,:)+90.0,xyzlatlon(2,:)); % *SUBFUNCTION*

% now make vectors xyzE pointing east at each point
% (compute at lat=equator for easier calculation)
[xyzE] = vsetV(zeros(1,size(xyzlatlon,2)),xyzlatlon(2,:)+90.0); % *SUBFUNCTION* 

% now assemble global transformation matrix RR that transforms all local
% coordinate directions (r,E,N) into global ones (x,y,z).
RRi  = repmat(uint32(1:ndof),3,1);
RRj  = [reshape(uint32(1:ndof)',3,[]);
        reshape(uint32(1:ndof)',3,[]);
        reshape(uint32(1:ndof)',3,[])];
RRv  = [rhat(1,:); rhat(2,:); rhat(3,:);
        xyzE(1,:); xyzE(2,:); xyzE(3,:);
        xyzN(1,:); xyzN(2,:); xyzN(3,:)];
RR   = sparse2(RRi(:),RRj(:),RRv(:)); % assembles sparse matrix from i,j,v triplets

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
  
  ind      = aa > 0.0;
  if any(ind)
      alat(ind)  = atan ( a(3,ind) ./ aa(ind)  ) .* rad2deg;
      along(ind) = atan2( a(2,ind), a(1,ind) ) .* rad2deg;
  end
  
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