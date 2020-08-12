function dofs = nod2dof(nods,dim)

[i,j] = size(nods);
switch dim
    case 1
        dofs = nods;
        return
    case 2
        dofs = [2.*nods(:)'-1; 2.*nods(:)'];
    case 3
        dofs = [3.*nods(:)'-2; 3.*nods(:)'-1; 3.*nods(:)'];
    otherwise
        error(' dim must be 1, 2 or 3');
end

dofs  = dofs(:);
if i==1
    dofs = dofs(:)';
elseif j==1
    dofs = dofs(:);
else
%     disp('Will return dofs as COLUMN vector!');
    dofs = dofs(:);
end

end % END OF FUNCTION nod2dof