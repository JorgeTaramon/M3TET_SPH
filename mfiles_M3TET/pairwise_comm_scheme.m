function comm_scheme = pairwise_comm_scheme(allSDNB)

verbose        = 0;
nsd            = size(allSDNB,2);
min_comm_steps = 1e8;
nNB_per_SD     = sum(allSDNB>0,1);

for iii=1:4
    switch iii
        case 1
            SDorder = 1:nsd;
            str = ' Order=1:nsd';
        case 2
            SDorder = nsd:-1:1;
            str = ' Order=nsd:-1:1';
        case 3
            [tmp,SDorder] = sort(nNB_per_SD,'descend');
            str = ' Order=nNB descending';
        case 4
            [tmp,SDorder] = sort(nNB_per_SD,'ascend');
            str = ' Order=nNB ascending';
    end

    comm_scheme   = zeros(2*nsd,nsd);
    ind           = zeros(nsd);
    for jsd=1:nsd
        isd = SDorder(jsd);
        nNB = length(find(allSDNB(:,isd)));
        for iNB=1:nNB
            NB     = allSDNB(iNB,isd);
            if ~ind(NB,isd)
%                 % This version is much slower than the 1-line command below
%                 my_ind  = find(comm_scheme(:,isd)==0);
%                 nb_ind  = find(comm_scheme(:,NB)==0);
%                 loc     = min(intersect(my_ind,nb_ind));
                % This version is much faster as it avoids repeated
                % "intersect" commands
                loc     = find( sum(comm_scheme(:,[isd NB]),2)==0,1,'first');
                comm_scheme(loc,isd) = NB;
                comm_scheme(loc,NB ) = isd;
                ind(isd,NB)          = 1;
            end
        end
    end
    max_nNB = find(sum(comm_scheme,2)>0,1,'last');
    
    if verbose
        disp(' Communication scheme: Subdomain...')
        disp(sprintf(' %2i',1:nsd));
        disp(' talks to...')
        for iNB=1:max_nNB
            disp(sprintf(' %2i',comm_scheme(iNB,:)));
        end
        disp(sprintf(' %s  -->  max_nNB=%2i',str,max_nNB));
        disp(' ');
    end

    % Check for minimum number of communication cycles
    if max_nNB<min_comm_steps
        min_comm_steps   = max_nNB;
        best_comm_scheme = comm_scheme(1:max_nNB,:);
    end
end
comm_scheme = best_comm_scheme;

end % END OF SUBFUNCTION find_pairwise_comm_scheme