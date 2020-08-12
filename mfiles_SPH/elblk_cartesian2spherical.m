function [theta,phi,r] = elblk_cartesian2spherical(x,y,z)
% Usage: [theta,phi,r] = elblk_cartesian2spherical(x,y,z)
%
% Purpose:
%   Transform Cartesian coordinates into spherical coordinates
%   
% Input:
%   x     : [matrix] : Cartesian x-coordinates (in km)
%   y     : [matrix] : Cartesian y-coordinates (in km)
%   z     : [matrix] : Cartesian z-coordinates (in km)

%
% Output:
%   theta : [matrix] : Theta of spherical coordinates (in rad)
%   phi   : [matrix] : Phi of spherical coordinates (in rad)
%   r     : [matrix] : Radius of spherical coordinates (in km)
%
% Colatitude (theta) is measured from +z axis to -z (range 0 to pi).
% Longitude (phi) is measured from +X axis in counterclowise direction
% (west to east) (range 0 to 2*pi).
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
% JMT Jun 2016

r            = sqrt(x.^2 + y.^2 + z.^2);
theta        = atan2(sqrt(x.^2 + y.^2),z);
phi          = atan2(y,x);
phi(phi < 0) = phi(phi < 0) + 2*pi;

end % END OF FUNCTION elblk_cartesian2spherical