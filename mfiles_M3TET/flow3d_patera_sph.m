function [VAR,Uinfo] = flow3d_patera_sph(VAR,MESH,COMM,UBC,SETTINGS,PHYSICS,NUMSCALE,dt)

if SETTINGS.nmg==1 && COMM.nsd==1
    [VAR,Uinfo] = flow3d_patera_nomg_sph(VAR,MESH,UBC,SETTINGS,PHYSICS,NUMSCALE,dt);
else
    [VAR,Uinfo] = flow3d_patera_mg_sph_p(VAR,MESH,COMM,UBC,SETTINGS,PHYSICS,NUMSCALE,dt);
end

end