function DensEl = calc_element_dens(VAR,SETTINGS,PHYSICS,EL2NOD,PhaseID,els,Nip,th_ip,ph_ip,r_ip)

nnodel = length(Nip);
switch SETTINGS.method_eval_dens
    case 'interp_nodal'
        % interpolation from nodal values to integration points
        DensEl = VAR.Dens( EL2NOD(1:nnodel,els) )'*Nip;

    case 'elem_phases'
        % use the constant densities defined in PHYSICS.Dens
        DensEl = PHYSICS.Dens(PhaseID(els))';

    case 'elem_var'
        DensEl = VAR.Dens(els);
        
    case 'fct_handle'
        DensEl = PHYSICS.Dens_fct(r_ip,th_ip);
        
    case 'custom01'
        % use the constant densities defined in PHYSICS.Dens and add
        % temperature effect
        DensEl = PHYSICS.Dens(PhaseID(els))';
        dT     = (VAR.T( EL2NOD(1:nnodel,els) )'*Nip) - PHYSICS.T_ref;
        DensEl = DensEl .* (1 - PHYSICS.alpha(PhaseID(els)) * dT);

    otherwise
        error(' Unknown "case" for evaluating density at integration points.');
end

% Make sure a column-vector is returned
DensEl = DensEl(:);

% DENSITY CUT-OFFS (IF DEFINED)
if isfield(PHYSICS,'minDens') && ~isempty(PHYSICS.minDens)
    DensEl = max(DensEl,PHYSICS.minDens); % lower cut-off for density
end
if isfield(PHYSICS,'maxDens') && ~isempty(PHYSICS.maxDens)
    DensEl = min(DensEl,PHYSICS.maxDens); % upper cut-off for density
end

end % END FUNCTION calc_element_dens