function WS = calc_tsearch2_ws(GCOORD,EL2NOD)

% Prepare element neighbor information for MUTILS' tsearch2. Providing this
% information for tsearch2 speeds up its performance. When running on older
% Matlab releases tsearch2 cannot calculate the neighbor information. In
% this case ("catch") it has to be calculated before calling tsearch2.

nVnod  = max(max(EL2NOD(1:4,:)));
xyz_pt = mean(GCOORD(:,EL2NOD(:,1)),2);
try
    opts.verbosity = 0;
    [~,WS]         = tsearch2(GCOORD(:,1:nVnod),EL2NOD([1 2 4 3],:),xyz_pt,[],[],opts);

catch
    try
        % The above command failed. In older Matlab releases tsearch2 
        % requires at least the element neighbor information. So I 
        % calculate WS.NEIGHBORS and try again.
        WS.NEIGHBORS = calc_tetra_neighbors(EL2NOD);
        % tsearch2 uses different local node numbering for tetrahedra:
        WS.NEIGHBORS = WS.NEIGHBORS([1 2 4 3],:);
        [~,WS] = tsearch2(GCOORD(:,nVnod),EL2NOD([1 2 4 3],:),xyzPT,WS);
    catch
        error('Calculating WS for tesearch2 fails');
    end
end

end