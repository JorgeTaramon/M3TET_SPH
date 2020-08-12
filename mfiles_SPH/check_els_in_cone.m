function [els_in_cone,els_within_cone,els_cone_bnd,els_outside_cone,els_in_cone_isoparametric] = ...
    check_els_in_cone(GCOORD_SPH,EL2NOD,theta_cone)
% Usage: [els_in_cone,els_within_cone,els_cone_bnd,els_outside_cone,els_in_cone_isoparametric] = ...
%   check_els_in_cone(GCOORD_SPH,EL2NOD,theta_cone)
% 
% Purpose: 
%   Check elements in relation to the cone
%
% Input:
%   GCOORD_SPH         : [matrix] : spherical coordinates for 4-node or 
%                                   10-node tetrahedrons mesh 
%   EL2NOD             : [matrix] : finite element connectivity matrix 
%                                   (4 x nel) (10 x nel)
%   theta_cone         : [scalar] : half aperture of the cone. 
%                                   0° < theta_cone < 90°
% 
% Output:
%   els_in_cone               : [vector] : elements in the cone (elements
%                                          within the cone + elements
%                                          crossing the cone boundary)
%   els_within_cone           : [vector] : elements within the cone
%   els_cone_bnd              : [vector] : elements crossing the cone
%                                          boundary 
%   els_outside_cone          : [vector] : elements outside the cone
%   els_in_cone_isoparametric : [vector] : elements crossing the cone
%                                          boundary and having at least 1
%                                          edge (bar) outside the cone
%
%                   FRONT VIEW                                            TOP VIEW                      
%
%                       Z                                                                             
%                       |                                                    |                         
%              _________|_________                                           |                         
%              \        |        /                                           |                         
%               \       |       /                                            |                         
%                \      |      /                                             |                         
%                 \2*theta_north                                             |                         
%                  \  __|__  /                                               |                         
%                   \/  |  \/                                             ___|___                         
%                    \  |  /                                             /   |   \                      
%                     \ | /                                             /    |    \                     
%                      \|/                                             |     |     |                    
%   --------------------|----------------------> Y       --------------------|----------------------> Y
%                      /|\                                             |     |     |                    
%                     / | \                                             \    |    /                     
%                    /  |  \                                             \___|___/                         
%                   /   |   \                                                |                         
%                  /    |    \                                               |                         
%                 /     |     \                                              |                         
%                /      |      \                                             |                         
%               /       |       \                                            |                         
%              /________|________\                                           |                         
%                       |                                                    |                         
%                                                                                                      
%                                                                            X                         
%
% JMT Jun 2016

theta_cone    = theta_cone*pi/180; % theta cone in rad
theta_north   = theta_cone;        % theta angle from the +Z axis to the generatrix
theta_south   = pi-theta_cone;     % theta angle from the -Z axis to the generatrix
nel           = size(EL2NOD,2);    % number of elements
theta         = GCOORD_SPH(1,1:max(max(EL2NOD(1:4,:))));
if nel > 1
    TH2EL = theta(EL2NOD(1:4,:));
else
    TH2EL = theta(EL2NOD(1:4,:))';
end

theta_inside_cone         = sum(TH2EL <= theta_north | TH2EL >= theta_south,1); % boolean vector for elements in the cone
els_in_cone               = find(theta_inside_cone); % indices for those elements in the cone boundary 
                                                     % (elements within the cone boundary + elements on the cone boundary)
els_within_cone           = find(theta_inside_cone == 4); % indices for elements within the cone boundary
els_outside_cone          = find(theta_inside_cone == 0); % indices for elements outside the cone boundary
els_cone_bnd              = find(theta_inside_cone > 0 & theta_inside_cone < 4); % indices for those elements crossing the cone boundary
els_in_cone_isoparametric = find(theta_inside_cone > 0 & theta_inside_cone < 3); % indices for those elements crossing the cone boundary and
                                                                                 % having at least 1 edge (bar) outside the boundary cone. 
                                                                                 % That means these elements will have spherical curved edges
                                                                                 % in the spherical rotated frame, so we need to track them 
                                                                                 % in order to use isoparametric spherical-to-local routine 
                                                                                 % when computing stiffness matrices
end % END OF FUNCTION check_els_in_cone