function [els_cross_2pi,VCOORD_ph] = check_phi(GCOORD_SPH,EL2NOD)
% Usage: [els_cross_2pi,VCOORD_ph] = check_phi(GCOORD_SPH,EL2NOD)
% 
% Purpose: Check which elements cross 2pi.
%
% Input:
%   GCOORD_SPH    : [matrix] : spherical coordinates for 4-node or 10-node 
%                              tetrahedrons mesh 
%   EL2NOD        : [matrix] : finite element connectivity matrix 
%                              (4 x nel) (10 x nel)
% 
% Output:
%   els_cross_2pi : [vector] : indices for those elements crossing 2pi
%   VCOORD_ph     : [matrix] : phi coordinates of the vertices for all
%                              elements after changing, for I elements, 
%                              those values of phi ~0 by ~0 + 2pi
%
% Check if all nodes within the same element are separated a distance less
% than pi --> i.e. check if there is any element whose nodes can be in the
% 1st and 5th quadrants (0 < phi < 90) and in the 4th and 8th quadrants 
% (270 < phi < 360). These nodes are close in Cartesian coordinates,
% however, they are ~2pi (360 degrees) far away each other in spherical
% coordinates.
% NOTE: phi (longitude) is measured from +X axis in counterclockwise 
% direction (range 0 to 2pi).
%
%               TOP VIEW
%
%                   |
%                   |
%                   |
%                   |
%                   |
%                   |
%                   |
%                   |
%   ----------------|----------------------> Y
%                   |
%                   |
%        1 *_- _ -_-|- - - _-* 3
%             - _  -|- * -  /
%                 - | 4 \  /
%                   |- _ \/
%                   |    * 2
%                   |
%                   |
%                   |
%                    
%                   X
%
% JMT Jun 2016

% % check if input GCOORD_SPH is actually in spherical coordinates
% if max(max(GCOORD_SPH(1,:)))/max(max(GCOORD_SPH(3,:))) > 0.01
%     % input GCOORD_SPH is in Cartesian coordiantes -> change to spherical coordinates
%     GCOORD_SPH = cartesian2spherical(GCOORD_SPH);
% end

nel                     = size(EL2NOD,2);          % number of elements
VCOORD_ph               = reshape(GCOORD_SPH(2,EL2NOD(1:4,1:nel)),4,nel); % phi coordiante of the vertices
ang_dist                = zeros(6,nel); 
ang_dist(1,:)           = abs(VCOORD_ph(1,:)-VCOORD_ph(2,:)); % angular distance between vertex 1 and vertex 2
ang_dist(2,:)           = abs(VCOORD_ph(2,:)-VCOORD_ph(3,:)); % angular distance between vertex 2 and vertex 3
ang_dist(3,:)           = abs(VCOORD_ph(3,:)-VCOORD_ph(1,:)); % angular distance between vertex 3 and vertex 1
ang_dist(4,:)           = abs(VCOORD_ph(1,:)-VCOORD_ph(4,:)); % angular distance between vertex 1 and vertex 4
ang_dist(5,:)           = abs(VCOORD_ph(2,:)-VCOORD_ph(4,:)); % angular distance between vertex 2 and vertex 4
ang_dist(6,:)           = abs(VCOORD_ph(3,:)-VCOORD_ph(4,:)); % angular distance between vertex 3 and vertex 4
ang_dist_longer_than_pi = ang_dist > pi; % angular distance longer than pi 
                                         % (means that some nodes are in the 1st/5th quadrant and some others nodes are in the 4th/8th quadrant)
els_cross_2pi           = sum(ang_dist_longer_than_pi,1); % boolean vector for elements crossing phi = 2pi
els_cross_2pi           = find(els_cross_2pi);            % indices of those elements crossing phi = 2pi

for i=1:size(els_cross_2pi,2)
    VCOORD_ph(VCOORD_ph(:,els_cross_2pi(i)) < pi,els_cross_2pi(i)) = VCOORD_ph(VCOORD_ph(:,els_cross_2pi(i)) < pi,els_cross_2pi(i)) + 2*pi;
    % change those values of phi ~0 by ~0 + 2pi
end

% if size(els_cross_2pi,2) <= 500
%     figure(60)
%     clf
%     GCOORD = spherical2cartesian(GCOORD_SPH);
%     lightGrey = 0.90*[1 1 1];
%     [x_sphere,y_sphere,z_sphere] = sphere(50);
%     x_sphere = x_sphere*6371;
%     y_sphere = y_sphere*6371;
%     z_sphere = z_sphere*6371;
%     axis equal
%     surface(x_sphere,y_sphere,z_sphere,'FaceColor','none','EdgeColor',lightGrey)
%     hold on
%     faceColor = [0.6875 0.8750 0.8984];
%     tetramesh(EL2NOD(:,els_cross_2pi)',GCOORD','FaceColor',faceColor,'FaceAlpha',0.3)
%     scatter3(GCOORD(1,EL2NOD(1:4,els_cross_2pi)),...
%              GCOORD(2,EL2NOD(1:4,els_cross_2pi)),...
%              GCOORD(3,EL2NOD(1:4,els_cross_2pi)),'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
%     axis([-6371 6371 -6371 6371 -6371 6371])
%     view(142.5,30)
%     grid on
%     faces = [1 2 3 4];
%     verts = [6371 0 6371; 0 0 6371; 0 0 -6371; 6371 0 -6371];
%     patch('Faces',faces,'Vertices',verts,'FaceColor','g','FaceAlpha',0.1)
%     xlabel('x (km)')
%     ylabel('y (km)')
%     zlabel('z (km)')
%     title('elements crossing \phi = 2\pi')
% end
end % END OF FUNCTION check_phi