function [RR] = make_rotation_matrix(GCOORD)
% Usage: [RR] = make_rotation_matrix(GCOORD)
% 
% Purpose: 
%   Create rotation matrix [RR] to transform spherical coordinates
%   (theta,phi,r) (or local coordinates (Colatitude, Longitude, Radial
%   direction)) to Cartesian coordinates (x,y,x). This is needed to create
%   a shear-stress free BC at the CMB (where only vR=0 will be a prescribed
%   BC) and to impose the velocity BC from plate tectonic reconstructions
%   at the Earth's surface.
%   The rotation matrix is the result after applying one rotation around
%   the Z axis and another one around the Y' axis.
%
%                 Z
%                 |
%                 |_
%                 | \
%                 |  \
%                 |   \
%                 |    \
%                 |     +
%                 |    /.
%                 |th / .
%                 |  /r .
%                 |^/   .
%                 |/    .
%                _-_----.-------_----> Y
%              _-\_/\   .     _-
%            _-  phi \  .   _-
%          _-         \ . _-
%        _-____________\.-
%      _-                
%    X      
%
%   1st rotation (Counterlockwise rotation around Z axis an angle phi):
%
%   [x ]   [ cos(phi)  -sin(phi)   0  ][x']
%   |y | = | sin(phi)   cos(phi)   0  ||y'|
%   [z ]   [    0          0       1  ][z']
%
%   2nd rotation (Counterlockwise around Y' axis an angle theta):
%
%   [x']   [ cos(theta)  0  sin(theta)][x"]
%   |y'| = |     0       1      0     ||y"|
%   [z']   [-sin(theta)  0  cos(theta)][z"]
%
%   The rotation matrix [RR] to change from local coordinates
%   (x",y",z") to global coordinates (x,y,z) is given by:
%   
%   [x]   [cos(phi)cos(theta)  -sin(phi)   cos(phi)sin(theta)][x"]
%   |y| = |sin(phi)cos(theta)   cos(phi)   sin(phi)sin(theta)||y"|
%   [z]   [   -sin(theta)         0            cos(theta)    ][z"]
%   
%   where [x"; y"; z"] is actually [theta; phi; r], so, e.g., for velocity:
%   
%   [U_x]   [cos(phi)cos(theta)  -sin(phi)   cos(phi)sin(theta)][U_theta]
%   |U_y| = |sin(phi)cos(theta)   cos(phi)   sin(phi)sin(theta)||U_phi  |
%   [U_z]   [   -sin(theta)         0            cos(theta)    ][U_r    ]
%
% Input:
%   GCOORD : [matrix]    : Cartesian coordinates (x,y,z) of the mesh
%
% Output:
%   RR     : [sparsemat] : Rotation matrix
%
% JMT Jun 2016
%

GCOORD_SPH             = cartesian2spherical(GCOORD);
theta                  = GCOORD_SPH(1,:);
phi                    = GCOORD_SPH(2,:);
nnod                   = size(GCOORD,2);

% Pre-allocate memory: number of non-zero slots in [RR]
nnz_T                  = 9*nnod;
KKi_T                  = zeros(nnz_T,1); 
KKj_T                  = KKi_T; 
KK1_T                  = ones(nnz_T,1);

% Diagonal slots (this will use up the first 3*nnod slots)
KKi_T(1:3*nnod)        = (1:3*nnod)';
KKj_T(1:3*nnod)        = (1:3*nnod)';
KK1_T(1:3:3*nnod-2)    = cos(phi).*cos(theta); % element (1,1) in the rotation matrix [RR]
KK1_T(2:3:3*nnod-1)    = cos(phi);             % element (2,2) in the rotation matrix [RR]
KK1_T(3:3:3*nnod  )    = cos(theta);           % element (3,3) in the rotation matrix [RR]

% Off-diagonal slots (this will use up the last 6*nnod slots)
KKi_T(3*nnod+1:4*nnod) = (1:3:3*nnod-2);
KKj_T(3*nnod+1:4*nnod) = (2:3:3*nnod-1);
KK1_T(3*nnod+1:4*nnod) = -sin(phi);            % element (1,2) in the rotation matrix [RR]

KKi_T(4*nnod+1:5*nnod) = (1:3:3*nnod-2);
KKj_T(4*nnod+1:5*nnod) = (3:3:3*nnod  );
KK1_T(4*nnod+1:5*nnod) = cos(phi).*sin(theta); % element (1,3) in the rotation matrix [RR]

KKi_T(5*nnod+1:6*nnod) = (2:3:3*nnod-1);
KKj_T(5*nnod+1:6*nnod) = (1:3:3*nnod-2);
KK1_T(5*nnod+1:6*nnod) = sin(phi).*cos(theta); % element (2,1) in the rotation matrix [RR]

KKi_T(6*nnod+1:7*nnod) = (2:3:3*nnod-1);
KKj_T(6*nnod+1:7*nnod) = (3:3:3*nnod  );
KK1_T(6*nnod+1:7*nnod) = sin(phi).*sin(theta); % element (2,3) in the rotation matrix [RR]

KKi_T(7*nnod+1:8*nnod) = (3:3:3*nnod  );
KKj_T(7*nnod+1:8*nnod) = (1:3:3*nnod-2);
KK1_T(7*nnod+1:8*nnod) = -sin(theta);          % element (3,1) in the rotation matrix [RR]

KKi_T(8*nnod+1:9*nnod) = (3:3:3*nnod  );
KKj_T(8*nnod+1:9*nnod) = (2:3:3*nnod-1);
KK1_T(8*nnod+1:9*nnod) = 0;                    % element (3,2) in the rotation matrix [RR]

RR = sparse2(KKi_T(:),KKj_T(:),KK1_T(:));

end % END OF FUNCTION make_rotation_matrix