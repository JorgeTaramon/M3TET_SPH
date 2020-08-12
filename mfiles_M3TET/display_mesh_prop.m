function display_mesh_prop(MESH,COMM,SETTINGS,NUMSCALE)

fidl = SETTINGS.fid_log;

fprintf(fidl,'\n');
fprintf(fidl,' Domain size             :   %7.2f x %7.2f x %7.2f %s\n',...
    MESH.Lx,MESH.Ly,MESH.Lz,NUMSCALE.unit_L);
fprintf(fidl,' Domain edges min(x,y,z) :   %7.2f , %7.2f , %7.2f %s\n',...
    MESH.xmin,MESH.ymin,MESH.zmin,NUMSCALE.unit_L);
fprintf(fidl,' Domain edges max(x,y,z) :   %7.2f , %7.2f , %7.2f %s\n',...
    MESH.xmax,MESH.ymax,MESH.zmax,NUMSCALE.unit_L);
fprintf(fidl,' Smallest element        :   %9.2f %s^3 , %7.2f %s \n',...
    min(MESH.vol_el),NUMSCALE.unit_L,min(MESH.len_el),NUMSCALE.unit_L);
fprintf(fidl,' Largest element         :   %9.2f %s^3 , %7.2f %s \n',...
    max(MESH.vol_el),NUMSCALE.unit_L,max(MESH.len_el),NUMSCALE.unit_L);
fprintf(fidl,' Number of MG levels     :  %10i\n',MESH.nmg);
fprintf(fidl,' Number of rock types    :  %10i\n', length(unique(MESH.PhaseID)));

if COMM.nsd>1
    fprintf(fidl,' Number of nodes in subdomain on all multigrid levels\n   ');
    fprintf(fidl,' %10i',COMM.nnod_SD);
    fprintf(fidl,'\n Number of elements in subdomain on all multigrid levels\n   ');
    fprintf(fidl,' %10i',COMM.nel_SD);
    fprintf(fidl,'\n Total number of nodes in domain on all multigrid levels\n   ');
    fprintf(fidl,' %10i',COMM.nnod_D);
    fprintf(fidl,'\n Total number of elements in domain on all multigrid levels\n   ');
    fprintf(fidl,' %10i',COMM.nel_D);
else
    nnod_mg = zeros(1,MESH.nmg);
    nel_mg  = zeros(1,MESH.nmg);
    for img=1:MESH.nmg
        nnod_mg(img) = max(max(MESH.EL2NOD{img}));
        nel_mg(img)  = size(MESH.EL2NOD{img},2);
    end
    fprintf(fidl,' Total number of nodes in domain on all multigrid levels\n   ');
    fprintf(fidl,' %10i',nnod_mg);
    fprintf(fidl,'\n Total number of elements in domain on all multigrid levels\n   ');
    fprintf(fidl,' %10i',nel_mg);
end
fprintf(fidl,'\n\n');

end