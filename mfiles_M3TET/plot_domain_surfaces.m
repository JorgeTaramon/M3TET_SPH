function plot_domain_surfaces(FigNo,MESH,ind_plot,var,view_angle)

if numlabs>1
    return
end
if nargin<4
    error(' call "plot_domain_surfaces(FigNo,MESH,ind_plot,VAR.*)"');
end
if nargin<5 || length(view_angle)~=2
    view_angle = [-40 15];
end
clim = [min(var) max(var)];
if clim(1)==clim(2)
    clim = [clim(1)-1e-8 clim(2)+1e-8];
end

nVnod = max(max(MESH.EL2NOD{1}(1:MESH.nvertx,:)));
nel   = size(MESH.EL2NOD{1},2);

if length(var)==nel
    [EL2NOD_face,els_face] = domain_face_connectivity...
        (MESH.EL2NOD{1}(1:MESH.nvertx,:),MESH.PointID,MESH.DB_indices,MESH.GCOORD,3);
elseif length(var)==nVnod
    EL2NOD_face = domain_face_connectivity...
        (MESH.EL2NOD{1}(1:MESH.nvertx,:),MESH.PointID,MESH.DB_indices,MESH.GCOORD,3);
else
    EL2NOD_face = domain_face_connectivity...
        (MESH.EL2NOD{1}(1:MESH.nvertx,:),MESH.PointID,MESH.DB_indices,MESH.GCOORD,3);
end
% EL2NOD_face is a 6-element cell that stores the surface connectivities of the domain
% EL2NOD_face{1} is domain bottom face,...
% EL2NOD_face{2} is domain front face,...
% EL2NOD_face{3} is domain right face,...
% EL2NOD_face{4} is domain back face,...
% EL2NOD_face{5} is domain left face,...
% EL2NOD_face{6} is domain top face

if ishandle(FigNo)
    axes(FigNo);
elseif isnumeric(FigNo)
    figure(FigNo);clf
else
    error('FigNo must be an integer >0 or an axis handle');
end
meshcol    = 'none';
for i=1:2
    if ~ismember(i,ind_plot)
        continue
    end
    EL2NOD_2D = EL2NOD_face{i};
    nel_2d    = size(EL2NOD_2D,2);
    nnodel_2d = size(EL2NOD_2D,1);
    el2V      = reshape((1:nnodel_2d*nel_2d)',nnodel_2d,nel_2d)';
    if isempty(var)
        trimesh(el2V,MESH.GCOORD(1,EL2NOD_2D),MESH.GCOORD(2,EL2NOD_2D),MESH.GCOORD(3,EL2NOD_2D),...
        'FaceColor','w','EdgeColor','k','FaceAlpha',1,'EdgeAlpha',0.2);
    else
        if length(var)==nel
            Vdata   = var(els_face{i});
            Vdata   = repmat(Vdata(:),1,3)';
            facecol = 'flat';
        else
            Vdata   = var(EL2NOD_2D);
            facecol = 'interp';
        end
        trimesh(el2V,MESH.GCOORD(1,EL2NOD_2D),MESH.GCOORD(2,EL2NOD_2D),MESH.GCOORD(3,EL2NOD_2D),...
            Vdata,'FaceColor',facecol,'EdgeColor',meshcol,'EdgeAlpha',0.2); %,'FaceAlpha',0.8,'EdgeAlpha',0.3);
    end
    hold on
%     xmin = min(min(MESH.GCOORD(1,EL2NOD_2D)));
%     xmax = max(max(MESH.GCOORD(1,EL2NOD_2D)));
%     ymin = min(min(MESH.GCOORD(2,EL2NOD_2D)));
%     ymax = max(max(MESH.GCOORD(2,EL2NOD_2D)));
%     zmin = min(min(MESH.GCOORD(3,EL2NOD_2D)));
%     zmax = max(max(MESH.GCOORD(3,EL2NOD_2D)));
%     if ismember(i,[2 4])
%         patch([xmin xmax xmax xmin]',[ymin ymin ymax ymax]',...
%               [zmin zmin zmax zmax]',[1 1 1 1]',...
%               'Facecolor','none','Edgecolor','k','Linewidth',1.0);
%     else
%         patch([xmin xmax xmax xmin]',[ymin ymin ymax ymax]',...
%               [zmin zmax zmax zmin]',[1 1 1 1]',...
%               'Facecolor','none','Edgecolor','k','Linewidth',1.0);
%     end
end

set(gca,'FontSize',8);
if ~isempty(clim)
    set(gca,'CLim',clim);
end
xlabel('x');ylabel('y');zlabel('z');
view(view_angle);
axis equal tight
% lighting phong
% lighting gouraud
if ~isempty(var)
    colormap(jet(100));
    h = colorbar;
    set(h,'Location','SouthOutside');
end
drawnow

end