function [EL2NOD_p1,PhaseID] = tetmesh_p2_to_p1(GCOORD,EL2NOD_p2,PhaseID)
% Usage: [EL2NOD_p1,PhaseID] = tetmesh_p2_to_p1(GCOORD,EL2NOD_p2,PhaseID)
%
% Purpose: Creates a linear (4-node) connectivity matrix from a quadratic
%          (10-node) connectivity matrix for a 3D tetrahedra mesh.
%
% Input:
%   EL2NOD_p2 : [matrix] : finite element connectivity matrix (10 x nel)
%   PhaseID   : [vector] : Phase-ID for each element (1 x nel)
%
% Output:
%   EL2NOD_p1 : [matrix] : finite element connectivity matrix (4 x 8*nel)
%   PhaseID   : [vector] : Phase-ID for each element (1 x 8*nel)
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Dec 2016
%

show_splits = 0;

nel      = size(EL2NOD_p2,2);
if isempty(GCOORD)
    isplit = 3*ones(1,nel);
else
    
    diag2nod = [5 7; 6 8; 9 10];
    L_diag   = zeros(3,nel);
    for i=1:3
        dx_diag     = GCOORD(1,EL2NOD_p2(diag2nod(i,1),:)) - GCOORD(1,EL2NOD_p2(diag2nod(i,2),:));
        dy_diag     = GCOORD(2,EL2NOD_p2(diag2nod(i,1),:)) - GCOORD(2,EL2NOD_p2(diag2nod(i,2),:));
        dz_diag     = GCOORD(3,EL2NOD_p2(diag2nod(i,1),:)) - GCOORD(3,EL2NOD_p2(diag2nod(i,2),:));
        L_diag(i,:) = sqrt(dx_diag.^2 + dy_diag.^2 + dz_diag.^2);
    end
    [~,isplit] = min(L_diag,[],1);
end
% isplit(:) = 3; disp('Setting p2-to-p1-split to 3 !!!!');

el2child    = cell(1,3);
% Split for diagonal inside parent connects nodes 5 and 7:
el2child{1} = [ 1  5  9  8 ...  % child 1
                5  2  6 10 ...  % child 2
                9  6  3  7 ...  % child 3
                8 10  7  4 ...  % child 4
                10 8  7  5 ...  % child 5
                6  9  5  7 ...  % child 6
                8  7  5  9 ...  % child 7
                10 5  7  6 ];   % child 8

% Split for diagonal inside parent connects nodes 6 and 8:
el2child{2} = [ 1  5  9  8 ...  % child 1
                5  2  6 10 ...  % child 2
                9  6  3  7 ...  % child 3
                8 10  7  4 ...  % child 4
                5  9  8  6 ...  % child 5
                10 7  6  8 ...  % child 6
                9  8  6  7 ...  % child 7
                5  6  8 10 ];   % child 8

% Split for diagonal inside parent connects nodes 9 and 10 (standard split)
el2child{3} = [ 1  5  9  8 ...  % child 1
                5  2  6 10 ...  % child 2
                9  6  3  7 ...  % child 3
                8 10  7  4 ...  % child 4
                8  5  9 10 ...  % child 5
                7  6 10  9 ...  % child 6
                5  9 10  6 ...  % child 7
                8 10  9  7 ];   % child 8

if show_splits
    % FigNo = 91; % helper for split with diag 5-7
    % el    = 2;
    % figure(FigNo);clf
    % plot_tetra(FigNo,GCOORD,EL2NOD_p2([2 4 3 1 10 7 9 5 6 8],:),el,4,'k','--',2);
    % view(-60,-85); hold all

    % FigNo = 92; % helper for split with diag 6-8
    % el    = 2;
    % figure(FigNo);clf
    % plot_tetra(FigNo,GCOORD,EL2NOD_p2([1 4 2 3 8 10 6 9 5 7],:),el,4,'k','--',2);
    % view(-60,-85); hold all

    % FigNo = 93; % helper for split with diag 9-10 (standard split)
    % el    = 2;
    % figure(FigNo);clf
    % plot_tetra(FigNo,GCOORD,EL2NOD_p2,el,4,'k','--',2);
    % view(-60,-85); hold all
    
    el    = 20;
    for j=1:3
        FigNo     = 1011;
        EL2NOD_p1 = reshape(EL2NOD_p2(el2child{j},el),4,[]);
        colors = lines(4);
        for i=5:8
            figure(FigNo+i);clf
            plot_tetra(FigNo+i,GCOORD,EL2NOD_p2,el,0,'k','--',2);
            view(-60,-85); hold all
            plot_tetra(FigNo+i,GCOORD,EL2NOD_p1,i,1,colors(i-4,:),'--',1);
        end
    end
end

EL2NOD_p1 = zeros(4,8*nel,'uint32');
for i=1:3
    iel_p2 = find(isplit==i);
    if ~isempty(iel_p2)
        iel_p1 = repmat(8*(iel_p2-1),8,1)+repmat((1:8)',1,length(iel_p2));
        EL2NOD_p1(:,iel_p1(:)') = reshape(...
            EL2NOD_p2(el2child{i},iel_p2),4,[]);
    end
end

if nargin==3 && length(PhaseID)==nel
    PhaseID = ones(8,1)*double(PhaseID(:)');
    PhaseID = int32(PhaseID(:)');
end

end % END OF FUNCION tetmesh_p2_to_p1