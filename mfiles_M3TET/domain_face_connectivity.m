function [EL2NOD_face,els_face] = domain_face_connectivity(EL2NOD,PointID,DB_indices,GCOORD,nodesface)
% Usage: [EL2NOD_face,els_face] = domain_face_connectivity(EL2NOD,PointID,DB_indices,GCOORD,nodesface)
%
% Purpose: Creates a 2D connectivity matrices for all six faces of the
% domain. The returned connectivity matrix for triangular elements are
% the faces of tetrahedral elements connected to the domain faces.
%
% Input:
%   EL2NOD     : [matrix] : connectivity matrix (nnodel x nel)
%   PointID    : [vector] : domain boundary index for each node (1 x nnod)
%   DB_indices : [cell]   : defines which domain edges/lines belong to which face
%   GCOORD     : [matrix] : coordinates of all nodes in coarse mesh (3 x nnod)
%   nodesface  : [scalar] : nodes per face (3 or 6)
%
% Output:
%   EL2NOD_face : [cell]  : connectivity matrix for each face
%   els_face    : [cell]  : element numbers connected to each face
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Mar 2013
% JH Feb 2015 : also returns PhaseID
% JH Feb 2016 : returns element number instead of PhaseID
% JMT Mar 2017 : added option to compute quadratic faces 
%

if nargin<4
    FigNo  = 0;
    GCOORD = [];
else
    FigNo  = 0;
end

EL2NOD_face = cell(1,6);
els_face    = cell(1,6);

%    108---211---107    % Coordinate system to define xmin, xmax, ymin, ymax, zmin, zmax
%     /:         /|     %     y
%  212 :  306  210|     %   /
%   /  :       /  |     %  /
%105---:209---106 |     % 0------- x
%  |  208  304|  207    % |
%  |   :      |   |     % |
%  |305:      |303|     % |
% 205  : 302 206  |     %-z
%  | 104---203|---103
%  |  /       |  /
%  |204  301  |202
%  |/         |/
%101--- 201---102

if FigNo;figure(FigNo);clf;hold on;end

if nodesface == 3
    for i=1:2
        nods = find(ismember(PointID,DB_indices{i}));
        [EL2NOD_face{i},els_face{i}] = make_face_connectivity_linear(EL2NOD,nods,FigNo,GCOORD);
    end
elseif nodesface == 6
    for i=1:2
        nods = find(ismember(PointID,DB_indices{i}));
        [EL2NOD_face{i},els_face{i}] = make_face_connectivity_quadratic(EL2NOD,nods,FigNo,GCOORD);
    end
else
    error('nodesface has to be either 3 or 6')
end

end % END OF FUNCTION domain_face_connectivity

% #########################################################################
%                               SUBFUNCTIONS
% #########################################################################

function [EL2NOD_face,els_face] = make_face_connectivity_linear(EL2NOD,nods,FigNo,GCOORD)

nnodel      = size(EL2NOD,1);
nedgenods   = 3;
els_face    = find(sum(ismember(EL2NOD,nods),1)==nedgenods);
nel         = length(els_face);
tmp         = EL2NOD(:,els_face);
tmp         = tmp(:);
[p1,ind]    = ismember(tmp,nods);
EL2NOD_face = reshape(uint32(nods(ind(p1>0))),nedgenods,nel);
[face,~]    = ind2sub([nnodel,nel],find(p1==0));
EL2NOD_face(:,face==1) = EL2NOD_face([1 3 2],face==1);
EL2NOD_face(:,face==2) = EL2NOD_face([3 1 2],face==2);
EL2NOD_face(:,face==3) = EL2NOD_face([1 3 2],face==3);
EL2NOD_face(:,face==4) = EL2NOD_face([1 2 3],face==4);

% ALTERNATIVE FORMULATION 1
% nface      = 4;
% nedgenods  = 3;
% els        = find(sum(ismember(EL2NOD,nods),1)==nedgenods);
% tmp        = EL2NOD(:,els);
% ind        = find(ismember(tmp,nods));
% [lnod,~]   = ind2sub(size(EL2NOD),ind);
% lnod       = reshape(lnod,3,[]);
% elnod2face = [2 4 3;
%               4 1 3;
%               1 4 2;
%               1 2 3];
% EL2NOD_face = zeros(nedgenods,length(els),'uint32');
% for iface=1:nface
%     ind = all(ismember(lnod,elnod2face(iface,:)));
%     EL2NOD_face(:,ind) = EL2NOD(elnod2face(iface,:),els(ind));
% end

% ALTERNATIVE FORMULATION 2
% The above block replaces the next loop
% EL2NOD_face = zeros(nedgenods,nel,'int32');
% for ii=1:nel
%     iel            = els(ii); % global element number    
%     [nod_glb,p1,~] = intersect(EL2NOD(:,iel),nods);
%     face           = setdiff(1:4,p1);
%     switch face
%         case 1 % conterclockwise order on face: 2 4 3 --> 1 3 2
%             nod_glb = EL2NOD([2 4 3],iel);
%         case 2 % conterclockwise order on face: 4 1 3 --> 3 1 2
%             nod_glb = EL2NOD([4 1 3],iel);
%         case 3 % conterclockwise order on face: 1 4 2 --> 1 3 2
%             nod_glb = EL2NOD([1 4 2],iel);
%         case 4 % conterclockwise order on face: 1 2 3 --> 1 2 3
%             nod_glb = EL2NOD([1 2 3],iel);
%     end
%     EL2NOD_face(:,ii) = nod_glb;
% end

if FigNo
    figure(FigNo);
    trimesh(EL2NOD_face',GCOORD(1,:),GCOORD(2,:),GCOORD(3,:),'FaceColor','k');
    axis equal tight
    colormap(jet);
%     colorbar % can crash Matlab
    view(-35,16);
end

end % END OF SUBFUNCTION make_face_connectivity_linear

function [EL2NOD_face,els_face] = make_face_connectivity_quadratic(EL2NOD,nods,FigNo,GCOORD)

nnodel      = size(EL2NOD(1:4,:),1);
nedgenods   = 6;
els_face    = find(sum(ismember(EL2NOD,nods),1)==nedgenods);
nel         = length(els_face);
tmp         = EL2NOD(:,els_face);
tmp         = tmp(:);
[p1,ind]    = ismember(tmp,nods);
EL2NOD_face = reshape(uint32(nods(ind(p1>0))),nedgenods,nel);
tmp2        = EL2NOD(1:4,els_face);
tmp2        = tmp2(:);
p2          = ismember(tmp2,nods);
[face,~]    = ind2sub([nnodel,nel],find(p2==0));
EL2NOD_face(:,face==1) = EL2NOD_face([1 3 2 6 5 4],face==1);
EL2NOD_face(:,face==2) = EL2NOD_face([3 1 2 6 4 5],face==2);
EL2NOD_face(:,face==3) = EL2NOD_face([1 3 2 6 5 4],face==3);
EL2NOD_face(:,face==4) = EL2NOD_face([1 2 3 4 5 6],face==4);

if FigNo
    figure(FigNo);
    trimesh(EL2NOD_face',GCOORD(1,:),GCOORD(2,:),GCOORD(3,:),'FaceColor','k');
    axis equal tight
    colormap(jet);
%     colorbar % can crash Matlab
    view(-35,16);
end

end % END OF SUBFUNCTION make_face_connectivity_quadratic