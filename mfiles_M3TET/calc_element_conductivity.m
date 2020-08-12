function Cond = calc_element_conductivity(VAR,OPTS_T,PHYSICS,EL2NOD,PhaseID,els,Nip)

switch OPTS_T.method_eval_cond
    case 'interp_nodal'
        % interpolation from nodal values to integration points
        Cond = VAR.Cond( EL2NOD(:,els) )'*Nip;
        
    case 'elem_phases'
        % use the constant conductivity defined in PHYSICS.K
        if length(PhaseID)==1
            Cond = PHYSICS.K(PhaseID) * ones(length(els),1);
        else
            Cond = PHYSICS.K(PhaseID(els))';
        end
        
    case 'elem_var'
        if size(EL2NOD,2)==4*length(VAR.Cond)
            Cond = VAR.Cond(ceil(els./4));
        else
            Cond = VAR.Cond(els);
        end
        
    otherwise
        error(' Unknown "case" for evaluating conductivity at integration points.');
end

% Make sure a column-vector is returned
Cond  = Cond(:);

end % END FUNCTION calc_element_conductivity