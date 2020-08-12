function plot_var_sph(FigNo,MESH,var,els)

figure(FigNo);clf;
p   = MESH.GCOORD';
if nargin<4 || isempty(els)
    t = MESH.EL2NOD{1}';
else
    t = MESH.EL2NOD{1}(:,els)';
end
tri = surftri(p,t);
h   = trimesh(tri,p(:,1),p(:,2),p(:,3),var);
colormap(jet(100));
axis equal
set(h,'EdgeColor','k','FaceColor','interp');
% colorbar
cmin = min(var);
cmax = max(var);
dc   = 1e-4*(cmax-cmin);
if dc==0
    dc = 1;
end
caxis([cmin-dc cmax+dc]);
xlabel('x (km)');ylabel('y (km)');zlabel('z (km)');

end