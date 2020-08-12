function RhoCp_blk = calc_element_RhoCp(VAR,OPTS_T,PHYSICS,EL2NOD,PhaseID,els,Nip)

switch OPTS_T.method_eval_dens
    case 'interp_nodal'
        % interpolation from nodal values to integration points
        Dens_blk = VAR.Dens( EL2NOD(:,els) )'*Nip;
        
    case 'elem_phases'
        % use the constant densities defined in PHYSICS.Dens
        if length(PhaseID)==1
            Dens_blk = PHYSICS.Dens(PhaseID) * ones(length(els),1);
        else
            Dens_blk = PHYSICS.Dens(PhaseID(els))';
        end
        
    case 'elem_var'
        if size(EL2NOD,2)==4*length(VAR.Dens)
            Dens_blk = VAR.Dens(ceil(els./4));
        else
            Dens_blk = VAR.Dens(els);
        end
        
    otherwise
        error(' Unknown "case" for evaluating density at integration points.');
end

switch OPTS_T.method_eval_Cp
    case 'interp_nodal'
        % interpolation from nodal values to integration points
        Cp_blk = VAR.Cp( EL2NOD(:,els) )'*Nip;
        
    case 'elem_phases'
        % use the constant Cp defined in PHYSICS.Cp
        if length(PhaseID)==1
            Cp_blk = PHYSICS.Cp(PhaseID) * ones(length(els),1);
        else
            Cp_blk = PHYSICS.Cp(PhaseID(els))';
        end
        
    case 'elem_var'
        if size(EL2NOD,2)==4*length(VAR.Cp)
            Cp_blk = VAR.Cp(ceil(els./4));
        else
            Cp_blk = VAR.Cp(els);
        end
        
    otherwise
        error(' Unknown "case" for evaluating heat capacity at integration points.');
end

% Make sure a column-vector is returned
RhoCp_blk = Dens_blk(:) .* Cp_blk(:);

end % END FUNCTION calc_element_RhoCp