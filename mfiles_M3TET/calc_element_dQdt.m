function dQdt_blk = calc_element_dQdt(VAR,OPTS_T,PHYSICS,EL2NOD,PhaseID,els,Nip)

switch OPTS_T.method_eval_dQdt
    case 'zero_dQdt'
        dQdt_blk = [];
        return
        
    case 'interp_nodal'
        % interpolation from nodal values to integration points
        dQdt_blk = VAR.dQdt( EL2NOD(:,els) )'*Nip;
        
    case 'elem_phases'
        % use the constant dQdt defined in PHYSICS.dQdt
        if length(PhaseID)==1
            dQdt_blk = PHYSICS.dQdt(PhaseID) * ones(length(els),1);
        else
            dQdt_blk = PHYSICS.dQdt(PhaseID(els))';
        end
        
    case 'elem_var'
        if size(EL2NOD,2)==4*length(VAR.dQdt)
            dQdt_blk = VAR.dQdt(ceil(els./4));
        else
            dQdt_blk = VAR.dQdt(els);
        end
        
    otherwise
        error(' Unknown "case" for evaluating heat capacity at integration points.');
end

% Make sure a column-vector is returned
dQdt_blk = dQdt_blk(:);

end % END FUNCTION calc_element_dQdt