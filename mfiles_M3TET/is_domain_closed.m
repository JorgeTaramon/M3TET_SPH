function is_closed = is_domain_closed(PointID,iUfix,GCOORD,EL2NOD)
% Usage: is_closed = is_domain_closed(PointID,iUfix,GCOORD,EL2NOD)
%
% Purpose: Check if boundary conditions allow for no unconstrained flow
%          in or out of the domain. In this case the pressure is defined up
%          to an arbitrary constant (problem becomes indefinite). In the CG
%          solution we will then force this constant to be zero, which 
%          allows CG to solve this positive-indefinite problem.
% Input:
%
% Output:
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)

% JH Jan 2011
%

is_closed = 1;

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

% x-dofs on boundary
idofX = find(ismember(PointID,[303 305 205:208 202 204 210 212 101:108]));
idofX = 3*idofX-2;
% y-dofs on boundary
idofY = find(ismember(PointID,[302 304 205:208 201 203 209 211 101:108]));
idofY = 3*idofY-1;
% z-dofs on boundary
idofZ = find(ismember(PointID,[301 306 201:204 209:212 101:108]));
idofZ = 3*idofZ;

%showUBC(82,GCOORD,EL2NOD,iUfix)

idofX = setdiff(idofX,iUfix);
idofY = setdiff(idofY,iUfix);
idofZ = setdiff(idofZ,iUfix);

if ~isempty(idofX)
    fprintf(' Domain open for flow in x-direction\n');
    fprintf(' Some nodes with the following PointIDs are not constrained:\n');
    fprintf(' %1i',unique(PointID(ceil(idofX/3)))); fprintf('\n');
    is_closed = 0;
end
if ~isempty(idofY)
    fprintf(' Domain open for flow in y-direction\n');
    fprintf(' Some nodes with the following PointIDs are not constrained:\n');
    fprintf(' %1i',unique(PointID(ceil(idofY/3)))); fprintf('\n');
    is_closed = 0;
end
if ~isempty(idofZ)
    fprintf(' Domain open for flow in z-direction\n');
    fprintf(' Some nodes with the following PointIDs are not constrained:\n');
    fprintf(' %1i',unique(PointID(ceil(idofZ/3)))); fprintf('\n');
    is_closed = 0;
end

% Visualize constrained node boundary conditions (for checking BCs)
FigNo = 0; % figure number (0 for no check)
if FigNo
    figure(FigNo);clf;
    
    scatter3(-100,-100,0,200,...
        'Marker','o','LineWidth',2,...
        'MarkerFaceColor','k','MarkerEdgeColor','k'); hold on
    scatter3(-100,-100,0,100,...
        'Marker','o','LineWidth',2,...
        'MarkerFaceColor','b','MarkerEdgeColor','b');
    scatter3(-100,-100,0,40,...
        'Marker','o','LineWidth',2,...
        'MarkerFaceColor','r','MarkerEdgeColor','r');
    legend('x-dof fixed','y-dof fixed','z-dof fixed');
    
    tetramesh(EL2NOD(:,1:4),GCOORD',...
              'EdgeColor','g',...
              'FaceAlpha',0.05,'FaceColor','k',...
              'LineStyle','-','LineWidth',1);
    xlabel('X');ylabel('Y');zlabel('Z');
    axis([min(GCOORD(1,:)) max(GCOORD(1,:)) ...
          min(GCOORD(2,:)) max(GCOORD(2,:)) ...
          min(GCOORD(3,:)) max(GCOORD(3,:))]);
    hold all
    for ii=1:length(iUfix)
        nod = ceil(iUfix(ii)/3);
        dof = 3-(3*nod-iUfix(ii));
        switch dof
            case 1
                scatter3(GCOORD(1,nod),GCOORD(2,nod),GCOORD(3,nod),200,...
                    'Marker','o','LineWidth',2,...
                    'MarkerFaceColor','none','MarkerEdgeColor','k');
            case 2
                scatter3(GCOORD(1,nod),GCOORD(2,nod),GCOORD(3,nod),100,...
                    'Marker','o','LineWidth',2,...
                    'MarkerFaceColor','none','MarkerEdgeColor','b');
            case 3
                scatter3(GCOORD(1,nod),GCOORD(2,nod),GCOORD(3,nod),40,...
                    'Marker','o','LineWidth',2,...
                    'MarkerFaceColor','none','MarkerEdgeColor','r');
        end
        text(GCOORD(1,nod),GCOORD(2,nod),GCOORD(3,nod)-40,...
             num2str(nod),'Fontsize',14,'Fontweight','bold');
    end
end

end % END OF FUNCTION is_domain_closed