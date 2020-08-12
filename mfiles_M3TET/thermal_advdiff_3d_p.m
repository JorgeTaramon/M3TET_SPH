function [VAR,WS_SLM] = thermal_advdiff_3d_p(VAR,MESH,COMM,PHYSICS,SETTINGS,NUMSCALE,TBC,dt,WS_SLM)

% TODO: advection:
% push points near top/bottom boundaries into mantle
% (calculate mid-face radius --> minimum shift inside)
% when element is found, use true point coords to calculate lc
% include adiabate and SeLa-correction

FigNo = 0;
% =====================================================================
% Thermal diffusion over half the time step (dt/2)
% =====================================================================
switch SETTINGS.T_solver
    case 'Thermal'
        VAR.T = thermal3d_p(VAR.T,VAR,MESH,COMM,PHYSICS,SETTINGS,NUMSCALE,TBC,dt/2);
    case 'Diffusion'
        VAR.T = diffusion3d_p(VAR.T,PHYSICS.kappa,MESH,COMM,SETTINGS,TBC,dt/2);
end
if FigNo
    ind_plot = [1 2 5];
    plot_domain_surfaces(FigNo,MESH,ind_plot,VAR.T,[-200 15]);
    title('T (diff)'); % caxis([1200 1400])
end
% compare_serial_parallel_3d(MESH,COMM,VAR.T,'T1',0);

% =====================================================================
% Thermal advection over the full time step dt (semi-Lagrange method)
% =====================================================================
switch SETTINGS.edge_type
    case 'straight'
        [VAR.T,WS_SLM]  = advect3d_slm_p...
            (MESH,COMM,VAR.Ux,VAR.Uy,VAR.Uz,VAR.T,dt,SETTINGS.OPTS_SLM,WS_SLM);
    case 'curved'
        [VAR.T,WS_SLM]  = advect3d_slm_p_sph...
            (MESH,COMM,VAR.Ux,VAR.Uy,VAR.Uz,VAR.T,dt,SETTINGS.OPTS_SLM,WS_SLM);
end
if FigNo
    plot_domain_surfaces(FigNo+1,MESH,ind_plot,VAR.T);
    title('T (diff+adv)'); % caxis([1200 1400])
end
%  compare_serial_parallel_3d(MESH,COMM,VAR.T,'T2',0);

% =====================================================================
% Thermal diffusion over second half the time step (dt/2)
% =====================================================================
switch SETTINGS.T_solver
    case 'Thermal'
        VAR.T = thermal3d_p(VAR.T,VAR,MESH,COMM,PHYSICS,SETTINGS,NUMSCALE,TBC,dt/2);
    case 'Diffusion'
        VAR.T = diffusion3d_p(VAR.T,PHYSICS.kappa,MESH,COMM,SETTINGS,TBC,dt/2);
end
if FigNo
    plot_domain_surfaces(FigNo+2,MESH,ind_plot,VAR.T);
    title('T (diff+adv+diff)'); % caxis([1200 1400])
end

end % END OF FUNCTION thermal_advdiff_3d_p