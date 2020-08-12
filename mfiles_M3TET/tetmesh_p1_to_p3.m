function [GCOORD,EL2NOD,PointID] = tetmesh_p1_to_p3(GCOORD4,EL2NOD4,nnodel)
% Usage: [GCOORD,EL2NOD,PointID] = tetmesh_p1_to_p3(GCOORD4,EL2NOD4,nnodel)
%
% Purpose: Creates a cubic (20- or 21-node) tetrahedra mesh from a linear
%          (4-node) hexahedra mesh.
%
% Input:
%   GCOORD4 : [matrix]  : coordinates of all nodes (3 x nnod)
%   EL2NOD4 : [matrix]  : finite element connectivity matrix (nnodel x nel)
%   nnodel  : [integer] : number of nodes per element to be created (20 or 21)
%
% Output:
%   GCOORD  : [matrix] : coordinates of all nodes in new mesh
%   EL2NOD  : [matrix] : finite element connectivity matrix (20 x nel or 21 x nel)
%   PhaseID : [vector] : Phase-ID for each element
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Mar 2016
%

if ~ismember(nnodel,[20 21])
    error('nnodel must be 20 or 21');
end

FigNo = 0;

% size of mesh
if size(EL2NOD4,1)~=4
    error('Expecting a 4-node conenctivity matrix');
end
nel     = size(EL2NOD4,2);         % number of elements
nnod4   = double(max(EL2NOD4(:))); % number of nodes
GCOORD4 = GCOORD4(:,1:nnod4);

if FigNo
    figure(FigNo);clf;
    for iel=1:nel
        show_elements_3d(GCOORD4,EL2NOD4,iel,'k',1);hold on
    end
    for inod=1:size(GCOORD4,2)
        text(GCOORD4(1,inod),GCOORD4(2,inod),GCOORD4(3,inod),...
            ['    ' num2str(inod)],...
            'HorizontalAlignment','center',...
            'Color','r','Fontsize',20,'Fontweight','bold');
    end
    view(58,16);
end

% Create a pointer edges defined by their end-nodes
% (e.g. bar 1 has end-nodes [1 2], bar 2 has [2 3],...)
% First create a 20-node hexahedron
edges = [ 1  2 ... % 1st edge of parent
          2  3 ... % 2nd edge of parent
          3  4 ... % 3rd edge of parent
          4  1 ... % 4th edge of parent
          1  3 ... % 5th edge of parent
          4  2 ];  % 6th edge of parent

EDGE2NOD = reshape( EL2NOD4(edges',:),2,[] )';

% Find the edges that are shared by neighboring elements and return a unique
% list of edges.
[EDGE2NOD,~,ib] = unique_keep_order(EDGE2NOD); % *SUBFUNCTION*
nedge  = size(EDGE2NOD,1); % number of unique edges (i.e. number of new edge 
                           % nodes that have to be generated)

EL2EDGE = reshape(ib,6,nel); % element to bar connectivity after doubles 
                             % have been merged (removed)

% Coordinates of the nodes defining each bar's end points
xBarEnds = reshape(GCOORD4(1,EDGE2NOD'),2,nedge);
yBarEnds = reshape(GCOORD4(2,EDGE2NOD'),2,nedge);
zBarEnds = reshape(GCOORD4(3,EDGE2NOD'),2,nedge);

% Create two new nodes on each edge
LxBar      = xBarEnds(2,:)-xBarEnds(1,:);
xBarThirds = [xBarEnds(1,:)+(1/3)*LxBar
              xBarEnds(1,:)+(2/3)*LxBar];
xBarThirds = xBarThirds(:)';
LyBar      = yBarEnds(2,:)-yBarEnds(1,:);
yBarThirds = [yBarEnds(1,:)+(1/3)*LyBar
              yBarEnds(1,:)+(2/3)*LyBar];
yBarThirds = yBarThirds(:)';
LzBar      = zBarEnds(2,:)-zBarEnds(1,:);
zBarThirds = [zBarEnds(1,:)+(1/3)*LzBar
              zBarEnds(1,:)+(2/3)*LzBar];
zBarThirds = zBarThirds(:)';

% storage for cubic order connectivity matrix
EL2NOD         = zeros(16,nel,'uint32'); % must be 16, face and bubble nodes added later
EL2NOD(1: 4,:) = EL2NOD4;
EL2NOD(5:16,:) = nnod4 + reshape([2*EL2EDGE(:)'-1; 2*EL2EDGE(:)'],12,[]);

% storage for quadratic order (20-node) node coordinates
nnod16                   = nnod4 + 2*nedge;
GCOORD                   = zeros(3,nnod16);
GCOORD(:,1:nnod4)        = GCOORD4;
GCOORD(1,nnod4+1:nnod16) = xBarThirds;
GCOORD(2,nnod4+1:nnod16) = yBarThirds;
GCOORD(3,nnod4+1:nnod16) = zBarThirds;



% =========================================================================
% Now correct local node numbering for all edge nodes
% =========================================================================
nvertx = 4;
% local coordinates of all edge nodes (nodes 5:16)
r = [2/3 1/3  0   0   0   0  1/3 2/3 2/3 1/3  0   0 ];
s = [1/3 2/3 2/3 1/3  0   0   0   0   0   0  1/3 2/3];
t = [ 0   0  1/3 2/3 2/3 1/3  0   0  1/3 2/3  0   0 ];
N = sf_dsf_tet([r;s;t],nvertx,'matrix');

xn_vertx = reshape(GCOORD(1,EL2NOD(1:nvertx,:)),nvertx,[]);
xn_check = N' * xn_vertx;
xn_mesh  = reshape(GCOORD(1,EL2NOD(5:16,:)),12,[]);
diffx    = abs(xn_check-xn_mesh)>1e-12;
yn_vertx = reshape(GCOORD(2,EL2NOD(1:nvertx,:)),nvertx,[]);
yn_check = N' * yn_vertx;
yn_mesh  = reshape(GCOORD(2,EL2NOD(5:16,:)),12,[]);
diffy    = abs(yn_check-yn_mesh)>1e-12;
zn_vertx = reshape(GCOORD(3,EL2NOD(1:nvertx,:)),nvertx,[]);
zn_check = N' * zn_vertx;
zn_mesh  = reshape(GCOORD(3,EL2NOD(5:16,:)),12,[]);
diffz    = abs(zn_check-zn_mesh)>1e-12;
for inod=1:2:12
    iel_flip = find(diffx(inod,:) | diffy(inod,:) | diffz(inod,:));
    if ~isempty(iel_flip)
        EL2NOD([inod+4 inod+5],iel_flip) = EL2NOD([inod+5 inod+4],iel_flip);
    end
end
% =========================================================================

if FigNo
    figure(FigNo);
    for inod=nnod4+1:nnod16
        text(GCOORD(1,inod),GCOORD(2,inod),GCOORD(3,inod),num2str(inod),...
            'HorizontalAlignment','center',...
            'Color','b','Fontsize',16,'Fontweight','bold');
    end
end

% =========================================================================
% Add a node on each of the four faces of the tetrahedra
% =========================================================================
[GCOORD,EL2NOD] = tetmesh_add_facenods(GCOORD,EL2NOD);

if FigNo
    figure(FigNo);
    nnod20 = size(GCOORD,2);
    for inod=nnod16+1:nnod20
        text(GCOORD(1,inod),GCOORD(2,inod),GCOORD(3,inod),num2str(inod),...
            'HorizontalAlignment','center',...
            'Color','g','Fontsize',16,'Fontweight','bold');
    end
end

if nnodel==21
    % =====================================================================
    % Add a node on each in the center of the tetrahedra
    % =====================================================================
    [GCOORD,EL2NOD] = tetmesh_add_bubblenods(GCOORD,EL2NOD(1:20,:),1);
    
    if FigNo
        figure(FigNo);
        nnod = size(GCOORD,2);
        for inod=nnod20+1:nnod
            text(GCOORD(1,inod),GCOORD(2,inod),GCOORD(3,inod),num2str(inod),...
                'HorizontalAlignment','center',...
                'Color','m','Fontsize',16,'Fontweight','bold');
        end
    
        figure(FigNo+1);clf
        show_hexahedron(GCOORD,EL2NOD,1,'k',1);
        show_hexahedron(GCOORD,EL2NOD,2,'b',1);
        show_hexahedron(GCOORD,EL2NOD,3,'r',1);
        show_hexahedron(GCOORD,EL2NOD,4,'g',1);

        figure(FigNo+2);clf
        show_hexahedron(GCOORD,EL2NOD,1,'r',2);
        show_hexahedron(GCOORD,EL2NOD,2,'b',2);
        show_hexahedron(GCOORD,EL2NOD,3,'r',2);
        show_hexahedron(GCOORD,EL2NOD,4,'g',2);
    end
end

% PointID = pointID_box(GCOORD);
PointID = [];

end % END OF FUNCTION tetmesh_p1_to_p3