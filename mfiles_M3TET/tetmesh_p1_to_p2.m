function [GCOORD,EL2NOD,PointID] = tetmesh_p1_to_p2(GCOORD_c,EL2NOD_c,PointID_c,DB_indices)
% Usage: [GCOORD,EL2NOD,PointID] = tetmesh_p1_to_p2(GCOORD_c,EL2NOD_c,PointID_c,DB_indices)
%
% Purpose: Creates a quadratic order (10-node) from a linear (4-node)
%          tetrahedra mesh by adding edge nodes to all elements.
%
% Input: 
%   GCOORD_c   : [matrix] : coordinates of all nodes (3 x nnod)
%   EL2NOD_c   : [matrix] : connectivity matrix (4 x nel)
%   PointID_c  : [vector] : domain boundary index for each node
%   DB_indices : [matrix] : domain boundary indices for each face of the domain
%
% Output:
%   GCOORD     : [matrix] : coordinates of all nodes new mesh (3 x ?)
%   EL2NOD     : [matrix] : connectivity matrix (10 x nel)
%   PointID    : [vector] : domain boundary index for each node (1 x ?)
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH March 2013
%

% Convert a linear-order (4-node) tetrahedron mesh into a quadratic-order 
% (10-node) tetrahedron mesh. To do so, bars are created that connect 
% the vertex nodes between which a new node has to be created. Then a 
% unique list of these bars is calculated to avoid multiple nodes in the
% same place. Onc the unique bar-list is obtained, node coordinates and
% the new connectivity matrix are calculated.

FigNo = 0;

% size of coarse mesh
nel_c  = size(EL2NOD_c,2);   % number of elements
nnod_c = size(GCOORD_c,2);   % number of nodes

% Create a pointer bars defined by their end-nodes
% (e.g. bar 1 has end-nodes [1 5], bar 2 has [5 2],...)
bar2node = reshape( EL2NOD_c([ 1  2 ... % 1st edge of parent
                               2  3 ... % 2nd edge of parent
                               3  4 ... % 3rd edge of parent
                               4  1 ... % 4th edge of parent
                               1  3 ... % 5th edge of parent
                               4  2 ... % 6th edge of parent
                                    ]',:),2,[] )';

% Find the bars that are shared by neighboring elements and return a unique
% list of bars.
[bar2node,~,ib] = find_unique_bars(bar2node); % *SUBFUNCTION*
% [bar2node2,~,ib2] = unique(sort(bar2node,2),'rows','stable'); % *SUBFUNCTION*
% isequal(ib,ib2)

nbars  = size(bar2node,1); % number of unique bars (i.e. number of new edge 
                           % nodes that have to be generated)

EL2BAR = reshape(ib,6,nel_c); % element to bar connectivity after doubles 
                              % have been merged (removed)

% Coordinates of the nodes defining each bar's end points
xBarEnds = reshape(GCOORD_c(1,bar2node'),2,[]);
yBarEnds = reshape(GCOORD_c(2,bar2node'),2,[]);
zBarEnds = reshape(GCOORD_c(3,bar2node'),2,[]);

% Create new node at each bar mid-point
xBarMids = 0.5*sum(xBarEnds,1);
yBarMids = 0.5*sum(yBarEnds,1);
zBarMids = 0.5*sum(zBarEnds,1);

% storage for quadratic order (10-node) connectivity matrix
EL2NOD                  = zeros(10,nel_c,'uint32'); 
EL2NOD(1: 4,:)          = EL2NOD_c;
EL2NOD(5:10,:)          = nnod_c+EL2BAR; clear EL2BAR
% storage for quadratic order (10-node) node coordinates
nnod                    = nnod_c + nbars;
GCOORD                  = zeros(3,nnod);
GCOORD(:,1:nnod_c)      = GCOORD_c;
GCOORD(1,nnod_c+1:nnod) = xBarMids;
GCOORD(2,nnod_c+1:nnod) = yBarMids;
GCOORD(3,nnod_c+1:nnod) = zBarMids;

% patch_tetrahedron(GCOORD_c,EL2NOD_c,1,60);
% view(-72,32); axis([0 0.1 0 0.1 0.88 1])
% for i=1:8
%     patch_tetrahedron(GCOORD,EL2NOD,i,60+i);
%     view(-72,32); axis([0 0.1 0 0.1 0.88 1])
% end

if nargin==4
    % Calculate boundary index for the new nodes
    PointID = calc_PointIDs(GCOORD,EL2NOD,PointID_c,DB_indices,FigNo); % *SUBFUNCTION*
else
    PointID = [];
end

end % END OF FUNCTION tetmesh_p1_to_p2

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function [edges_uniq,IIu,JJu] = find_unique_bars(edges)
    if size(edges,2)~=2
        error('expecting a nx2 matrix');
    end

    % find FIRST occurrence of unique rows in edges
    % use sort so that vertices defining the edges are sorted increasingly
    [~,iu,ju]   = unique(sort(edges,2),'rows','first');

    % create index IIu to extract edges_uniq from edges
    [IIu,~,JJu] = unique(iu(ju));

    % extract unique edges from input
    edges_uniq  = edges(IIu,:);
end % END OF SUBFUNCTION find_unique_bars

% #########################################################################

function PointID = calc_PointIDs(GCOORD,EL2NOD,PointID_c,DB_indices,FigNo)

nVnod = max(max(EL2NOD(1: 4,:)));
nnod  = max(max(EL2NOD(5:10,:)));
nEnod = nnod-nVnod;

% Convert to old format (have to rewrite this function to use
% PointID-format)
% ===========================================================
DBnods_c      = find(PointID_c>0)';
if isempty(DBnods_c)
    PointID   = zeros(1,nnod,'int32');
    return
end
DBnods_c(:,2) = PointID_c( DBnods_c );

% Get boundary index information for new nodes
% ================================================
% 1) Create a pointer that for each edge node in the mesh returns the 2
%    bounding vertex nodes
pointer_Edge2Vnod                       = zeros(nEnod,2,'int32');
pointer_Edge2Vnod(EL2NOD( 5,:)-nVnod,:) = EL2NOD([1 2],:)';
pointer_Edge2Vnod(EL2NOD( 6,:)-nVnod,:) = EL2NOD([2 3],:)';
pointer_Edge2Vnod(EL2NOD( 7,:)-nVnod,:) = EL2NOD([3 4],:)';
pointer_Edge2Vnod(EL2NOD( 8,:)-nVnod,:) = EL2NOD([4 1],:)';
pointer_Edge2Vnod(EL2NOD( 9,:)-nVnod,:) = EL2NOD([1 3],:)';
pointer_Edge2Vnod(EL2NOD(10,:)-nVnod,:) = EL2NOD([2 4],:)';


% (2) Find fine mesh vertex nodes (i.e. coarse mesh nodes) that are in the
%     list of coarse domain boundary nodes (1st column in DBnods_c)
[tf,loc] = ismember(pointer_Edge2Vnod,DBnods_c(:,1)); % find DBnod_c in pointer_Edge2Vnod
tf       = find(tf);
% Create new pointer to the domain boundary index (2nd column in DBnods_c)
pointer_Edge2DBnod     = zeros(nEnod,2,'int32'); % create new pointer
pointer_Edge2DBnod(tf) = DBnods_c(loc(tf),2);

% (3)  Use the pointer pointer_Edge2DBnod to check (for each NEW node) 
%      what the domain boundary index of its two neighborung vertex nodes
%      are. The domain boundary index for the new nodes can be found using
%      the following logic:
%      Edge nodes cannot be generated in domain corners but only on
%      domain edges or faces.
%      Face nodes cannot be generated on domain edges (and of course not in
%      domain corners).
%
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

% Faces: 301 bottom
%        302 front
%        303 right
%        304 back
%        305 left
%        306 top
% Edges: 201 lower front
%        202 lower right
%        203 lower back
%        204 lower left
%        205 front left
%        206 front right
%        207 back right
%        208 back left
%        209 top front
%        210 top right
%        211 top back
%        212 top left
% Corners: 101 lower front left
%          102 lower front right
%          103 lower back right
%          104 lower back left
%          105 top front left
%          106 top front right
%          107 top back right
%          108 top back left
%
% special: 310 face nodes west of MAR
%          220 edge nodes west of MAR
%          311 face nodes east of MAR
%          221 edge nodes east of MAR
%          222 nodes on MAR
%          223 nodes in TF
%          110 nodes at JCT (intersection MAR-TF)

% There are 9 standard cases and 6 special cases in 3D:
%       1) both on same F ==> F
%       2) both on same E ==> E
%       3) C and E ==> E (C+E must be on same domain edge)
%       4) C and F ==> F (C+F must be on same domain face)
%       5) E and F ==> F (E+F must be on same domain face)
%       6) C and 0 ==> 0
%       7) E and 0 ==> 0
%       8) F and 0 ==> 0
%       9) 0 and 0 ==> 0 
%      S1) C1 and C2 ==> E between C1 and C2
%      S2) E1 and E2 ==> F between E1 and E2
%      S3) F1 and F2 ==> 0 (inside domain)
%      S4) C1 and E2 ==> 0 (inside domain)
%      S5) C1 and F2 ==> 0 (inside domain)
%      S6) E1 and F2 ==> 0 (inside domain)
%      The numbers in the special cases indicate DIFFERENT domain faces !!!

% % % For testing and understanding the logic:
% % pointer_Edge2DBnod = [301 301; % case 1
% %                        201 201; % case 2
% %                        101 201; % case 3
% %                        101 301; % case 4
% %                        201 301; % case 5
% %                        101  0 ; % case 6
% %                        201  0 ; % case 7
% %                        301  0 ; % case 8
% %                         0   0 ; % case 9
% %                        101 102; % special case 1
% %                        201 202; % special case 2
% %                        301 302];% special case 3
% %                        101 202; % special case 4
% %                        101 303; % special case 5
% %                        201 303];% special case 6

% Cases 1:9 are taken care of by chosing the LARGER index of the 2
% bounding vertex nodes. Treat all nodes as standard cases first:
DBindx  = max( pointer_Edge2DBnod,[],2 );
% HOWEVER, if any index is zero, the node will be inside the domain. Set
% these indices equal to zero now.
ind_0   = any(pointer_Edge2DBnod==0,2);
DBindx(ind_0) = 0;

% Need special treatment for the special cases...
% (a) all potential special cases
indS  = find( pointer_Edge2DBnod(:,1)~=pointer_Edge2DBnod(:,2) & ...
              pointer_Edge2DBnod(:,1)>0 & pointer_Edge2DBnod(:,2)>0 );
% (b) special case 1 (between 2 corners); all possible combinations
indS1 = indS(all(pointer_Edge2DBnod(indS,:)>100,2) & all(pointer_Edge2DBnod(indS,:)<200,2));
if ~isempty(indS1)
    ind_0          = indS1;
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[101 102]) ,2) );
    DBindx(indS1a) = 201; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[102 103]) ,2) );
    DBindx(indS1a) = 202; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[103 104]) ,2) );
    DBindx(indS1a) = 203; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[101 104]) ,2) );
    DBindx(indS1a) = 204; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[101 105]) ,2) );
    DBindx(indS1a) = 205; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[102 106]) ,2) );
    DBindx(indS1a) = 206; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[103 107]) ,2) );
    DBindx(indS1a) = 207; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[104 108]) ,2) );
    DBindx(indS1a) = 208; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[105 106]) ,2) );
    DBindx(indS1a) = 209; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[106 107]) ,2) );
    DBindx(indS1a) = 210; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[107 108]) ,2) );
    DBindx(indS1a) = 211; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[105 108]) ,2) );
    DBindx(indS1a) = 212; ind_0 = setdiff(ind_0,indS1a);
    
    for iDface=1:length(DB_indices)
        iSameFace = sum(ismember( pointer_Edge2DBnod(ind_0,:),DB_indices{iDface} ),2)==2;
        DBindx(ind_0(iSameFace)) = max(DB_indices{iDface});
        ind_0(iSameFace)         = [];
    end
    
    DBindx(ind_0)  = 0;
end

% (c) special case 2 (between 2 different edges); all possible combinations
indS2 = indS(all(pointer_Edge2DBnod(indS,:)>200,2) & all(pointer_Edge2DBnod(indS,:)<300,2));
if ~isempty(indS2)
    ind_0          = indS2;
    indS2a         = indS2( all( ismember(pointer_Edge2DBnod(indS2,:),[201 202 203 204]) ,2) );
    DBindx(indS2a) = 301; ind_0 = setdiff(ind_0,indS2a);
    indS2a         = indS2( all( ismember(pointer_Edge2DBnod(indS2,:),[201 205 206 209]) ,2) );
    DBindx(indS2a) = 302; ind_0 = setdiff(ind_0,indS2a);
    indS2a         = indS2( all( ismember(pointer_Edge2DBnod(indS2,:),[202 206 207 210]) ,2) );
    DBindx(indS2a) = 303; ind_0 = setdiff(ind_0,indS2a);
    indS2a         = indS2( all( ismember(pointer_Edge2DBnod(indS2,:),[203 207 208 211]) ,2) );
    DBindx(indS2a) = 304; ind_0 = setdiff(ind_0,indS2a);
    indS2a         = indS2( all( ismember(pointer_Edge2DBnod(indS2,:),[204 205 208 212]) ,2) );
    DBindx(indS2a) = 305; ind_0 = setdiff(ind_0,indS2a);
    indS2a         = indS2( all( ismember(pointer_Edge2DBnod(indS2,:),[209 210 211 212]) ,2) );
    DBindx(indS2a) = 306; ind_0 = setdiff(ind_0,indS2a);
    
    for iDface=1:length(DB_indices)
        iSameFace = sum(ismember( pointer_Edge2DBnod(ind_0,:),DB_indices{iDface} ),2)==2;
        DBindx(ind_0(iSameFace)) = max(DB_indices{iDface});
        ind_0(iSameFace)         = [];
    end
    
    DBindx(ind_0)  = 0;
end

% (d) special case 3
indS3 = indS(all(pointer_Edge2DBnod(indS,:)>300,2) & all(pointer_Edge2DBnod(indS,:)<400,2));

% (e) special cases 4, 5 and 6
indS456 = setdiff(indS,[indS1(:);indS2(:);indS3(:)]);
% Now loop over the domain faces and check if the different DB indices
% (e.g. C1 and E2) are on the same domain face as defined by the cell-matrix 
% "DB_indices". If they are on the same face, they are NOT INSIDE the
% domain!
for iDface=1:length(DB_indices)
    iSameFace = sum(ismember( pointer_Edge2DBnod(indS456,:),DB_indices{iDface} ),2)==2;
    indS456(iSameFace) = [];
end
% old version:
% Dfaces  = [301 201 202 203 204 101 102 103 104; % domain bottom     (zmin)
%            302 201 205 206 209 101 102 105 106; % domain front      (ymin)
%            303 202 206 207 210 102 103 106 107; % domain right side (xmax)
%            304 203 207 208 211 103 104 107 108; % domain back       (ymax)
%            305 204 205 208 212 101 104 105 108; % domain left side  (xmin)
%            306 209 210 211 212 105 106 107 108];% domain top        (zmax)
% % Now loop over the domain faces and check if the different DB indices
% % (e.g. C1 and E2) are on the same domain face as defined by matrix Dfaces.
% for iDface=1:size(Dfaces,1)
%     iSameFace = sum(ismember( pointer_Edge2DBnod(indS456,:),Dfaces(iDface,:) ),2)==2;
%     indS456(iSameFace) = [];
% end

% Special case 3,4,5 and 6 (between 2 different domain faces thus INSIDE
% the domain) ===> simply remove these nodes from the list of boundary nodes
inod_DB = (nVnod+1:nnod)';
DBindx ([indS3(:);indS456(:)]) = [];
inod_DB([indS3(:);indS456(:)]) = [];

% Remove all nodes that have a zero DB index
ind_0          = ~DBindx;
DBindx(ind_0)  = [];
inod_DB(ind_0) = [];

% merge coarse-mesh and new domain boundary nodes
DBnods                = [DBnods_c; [inod_DB DBindx]];
PointID              = zeros(1,nnod,'int32');
PointID(DBnods(:,1)) = DBnods(:,2);

% FigNo = 11;
if FigNo
    sfigure(FigNo);clf;
    
    uniqDBindx        = unique(DBnods(:,2));
    nDBindx           = length(uniqDBindx);
    col               = zeros(max(DBnods(:,2)),3);
    col(uniqDBindx,:) = lines(nDBindx);

    tetramesh(EL2NOD(1:4,:)',GCOORD',...
              'EdgeColor','g',...
              'FaceAlpha',0.6,'FaceColor','k',...
              'LineStyle','-','LineWidth',1);
    hold all
    xlabel('X');ylabel('Y');zlabel('Z');
    for ii=1:length(DBnods)
        nod  = DBnods(ii,1);
        indx = DBnods(ii,2);
        text(GCOORD(1,nod),GCOORD(2,nod),GCOORD(3,nod),num2str(indx),...
             'Fontweight','bold','Fontsize',14,'Color',col(indx,:));
        scatter3(GCOORD(1,nod),GCOORD(2,nod),GCOORD(3,nod),5,...
                 'Marker','o','LineWidth',2,...
                 'MarkerFaceColor',col(indx,:),'MarkerEdgeColor',col(indx,:));
    end
    
    sfigure(FigNo+1);clf;
    tetramesh(EL2NOD(1:4,:)',GCOORD',...
              'EdgeColor','g',...
              'FaceAlpha',0.6,'FaceColor','k',...
              'LineStyle','-','LineWidth',1);
    hold all
    xlabel('X');ylabel('Y');zlabel('Z');
    for ii=1:length(DBnods)
        nod  = DBnods(ii,1);
        indx = DBnods(ii,2);
        text(GCOORD(1,nod),GCOORD(2,nod),GCOORD(3,nod),num2str(nod),...
             'Fontweight','bold','Fontsize',14,'Color',col(indx,:));
        scatter3(GCOORD(1,nod),GCOORD(2,nod),GCOORD(3,nod),5,...
                 'Marker','o','LineWidth',2,...
                 'MarkerFaceColor',col(indx,:),'MarkerEdgeColor',col(indx,:));
    end
end

end % END OF SUBFUNCTION calc_PointIDs