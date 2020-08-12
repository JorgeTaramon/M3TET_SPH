function GCOORD = adjust_surface_radius(GCOORD,PointID,SurfID)

FigNo = 0;

nod_surf  = find(PointID==SurfID);
[az,el,r] = cart2sph(GCOORD(1,nod_surf),GCOORD(2,nod_surf),GCOORD(3,nod_surf));

if FigNo
    figure(FigNo); clf
    scatter3(GCOORD(1,nod_surf),GCOORD(2,nod_surf),GCOORD(3,nod_surf),10,r);
    colorbar
    axis equal tight
end

r0                 = max(r);
[x,y,z]            = sph2cart(az,el,r0);
GCOORD(1,nod_surf) = x;
GCOORD(2,nod_surf) = y;
GCOORD(3,nod_surf) = z;

if FigNo
    r = sqrt( x.^2 + y.^2 + z.^2 );
    figure(2); clf
    scatter3(x,y,z,20,r,'filled');
    colorbar
    axis equal tight
end

end