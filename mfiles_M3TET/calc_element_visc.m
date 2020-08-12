function ViscEl = calc_element_visc(VAR,SETTINGS,PHYSICS,EL2NOD,PhaseID,els,Nip,th_ip,ph_ip,r_ip)

nnodel = length(Nip);
switch SETTINGS.method_eval_visc
    case 'interp_nodal'
        % interpolation from nodal values to integration points
        ViscEl = VAR.Visc( EL2NOD(1:nnodel,els) )'*Nip;

    case 'mean_nodal'
        % harmonic mean for each element
        ViscEl = VAR.Visc( EL2NOD(1:nnodel,els) );
        ViscEl = nnodel./sum(1./ViscEl,1)';

    case 'elem_phases'
        % use the constant viscosities defined in PHYSICS.Visc
        ViscEl = PHYSICS.Visc( PhaseID(els) );

    case 'elem_var'
        ViscEl = VAR.Visc(els);
        
    case 'fct_handle'
        ViscEl = PHYSICS.Visc_fct(r_ip);

    case 'custom01'
        % PhaseID plus T-dependence
        ViscEl = PHYSICS.Visc( PhaseID(els) )';
        TempEl = VAR.T( EL2NOD(1:nnodel,els) )'*Nip;
        
        % Temperature effect on viscosity
        R      = 8.314472; % J/mol K; univ. gas constant
        T0_K   = PHYSICS.T_ref+273;
        eta    = PHYSICS.Ea(PhaseID(els))./(R*T0_K); 
        ViscEl = ViscEl .* exp( eta.*( (T0_K./(TempEl+273)) - 1) );
        
    otherwise
        error(' Unknown "case" for evaluating viscosity at integration points.');
end

% Make sure a column-vector is returned
ViscEl = ViscEl(:);

% VISCOSITY CUT-OFFS (IF DEFINED)
if isfield(PHYSICS,'minVisc') && ~isempty(PHYSICS.minVisc)
    ViscEl = max(ViscEl,PHYSICS.minVisc); % lower cut-off for viscosity
end
if isfield(PHYSICS,'maxVisc') && ~isempty(PHYSICS.maxVisc)
    ViscEl = min(ViscEl,PHYSICS.maxVisc); % upper cut-off for viscosity
end

% Check if that viscosities are within valid range
if any(ViscEl<=0)
    error(' Calculated a zero or negative viscosity. STOPPING.');
end

end % END OF FUNCTION calc_element_visc