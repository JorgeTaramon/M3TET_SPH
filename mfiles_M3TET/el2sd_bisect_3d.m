function [el2sd,topo] = el2sd_bisect_3d(GCOORD,EL2NOD,COMM,fewer_NB_wght,fid_log)
% Usage: [el2sd,topo] = el2sd_bisect_3d(GCOORD,EL2NOD,COMM,fewer_NB_wght,fid_log)
%
% Purpose: 
%
% Input:
%
% Output:
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH March 2013
%

% Get all possible bisections for given nsd
nsplit = domain_decomp(COMM.nsd); % *SUBFUNCTION*

% TRANSFORM CARTESIAN TO SPHERICAL COORDINATES (could be useful to find a 
% more efficient subdomain configuration in a sphere)
[azim,elev,rad] = cart2sph(GCOORD(1,:),GCOORD(2,:),GCOORD(3,:));
GCOORD          = [azim;elev;rad];

% Find best bisection by calculating the number of nodes shared between SDs
if numlabs==1
    % Serial version (all SDs test all splits)
    [el2sd,topo] = find_best_split(nsplit,GCOORD,EL2NOD',COMM,fewer_NB_wght,fid_log); % *SUBFUNCTION*
else
    % Parallel version (each SD tests different splits)
    [el2sd,topo] = find_best_split_p(nsplit,GCOORD,EL2NOD',COMM,fewer_NB_wght,fid_log); % *SUBFUNCTION*
end

end % END OF FUNCTION calculate_el2sd_bisect_3D

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function nsplit = domain_decomp(nsd)

[nx,ny,nz] = meshgrid(1:nsd,1:nsd,1:nsd);
ind        = nx.*ny.*nz == nsd;
nsplit     = [nx(ind) ny(ind) nz(ind)];
return

% The following does not give ALL possible configurations...
n = log2(nsd);
if n==round(n)
    iexp   = 1:n;
    k      = [1 2.^iexp];
    nsplit = [k(end:-1:1)' k' ones(n+1,1)];
    test   = prod(nsplit,2);
    if any(test~=nsd)
        error('domain decomposition failed');
    end
else
    nsplit(1,:) = [nsd 1 1];
    ind         = 1;
    ind1        = 1;
    while ind1<=ind
        nn   = nsplit(ind1,1);
        div1 = setdiff(unique(gcd(nn,1:nn)),[1 nn]);
        if ~isempty(div1)
            div1 = div1(1:ceil(length(div1)/2));
        end

        for id1=1:length(div1)
            test            = [nsplit(ind1,1)/div1(id1) nsplit(ind1,2)*div1(id1) nsplit(ind1,3)];
            nsplit(ind+1,:) = test;
            ind             = ind + 1;

            nn              = nsplit(ind,2);
            div2            = setdiff(unique(gcd(nn,1:nn)),[1 nn]);
            if ~isempty(div2)
                div2 = div2(1:ceil(length(div2)/2));
            end
            ind2 = ind;
            for id2=1:length(div2)
                nsplit(ind+1,:) = [nsplit(ind2,1) nsplit(ind2,2)/div2(id2) nsplit(ind2,3)*div2(id2)];
                ind             = ind + 1;
            end
        end
        ind1 = ind1 + 1;
    end

    nsplit = sort(nsplit,2,'descend');
    nsplit = unique(nsplit,'rows');

    test = nsplit(:,1) .* nsplit(:,2) .* nsplit(:,3);
    if any(test~=nsd)
        error('domain decomposition failed');
    end
end

m      = perms([1 2 3]);
nsplit = reshape(nsplit(:,m),[],3);
nsplit = unique(nsplit,'rows');

end % END OF SUBFUNCTION domain_decomp

% #########################################################################

function [el2sd,topo] = find_best_split(nsplit,GCOORD,EL2NOD,COMM,fewer_NB_wght,fid_log)

ntest = size(nsplit,1);
nnod  = double(max(max(EL2NOD)));
nsd   = COMM.nsd;

fprintf(fid_log,' Testing %1i configurations to find most efficient domain\n decomposition for %1i subdomains.\n',ntest,nsd);
fprintf(fid_log,' *********************************************************\n');
fprintf(fid_log,' case | # bi-sections |  # shared nodes  | # comm cycles \n');

pct_shared    = zeros(ntest,1);
ncomm         = zeros(ntest,1);
sumSDnods     = zeros(ntest,1);
invalid_split = zeros(ntest,1);
for itest=1:ntest
    nsplit_now = nsplit(itest,:);
    el2sd      = assign_el2sd(GCOORD,EL2NOD,nsd,nsplit_now); % *SUBFUNCTION*
    
    if length(unique(el2sd))<nsd
        invalid_split(itest)=1;
    end
    
    allSDNB = zeros(nsd-1,nsd);
    for isd=1:nsd
        NBs = el2sd( sum(ismember(EL2NOD(:,1:4),EL2NOD(el2sd==isd,1:4)),2)>0 );
        NBs = unique(NBs(NBs~=isd));
        nNB = length(NBs);
        allSDNB(1:nNB,isd) = NBs(:);
    end
    
    comm_scheme  = pairwise_comm_scheme(allSDNB);
    ncomm(itest) = size(comm_scheme,1);
     
    nnodSD = zeros(nsd,1);
    for isd=1:nsd
        nnodSD(isd) = length( unique( EL2NOD(el2sd==isd,: ) ) );
    end
    sumSDnods(itest)  = sum(nnodSD);
    pct_shared(itest) = 100*(sumSDnods(itest)-nnod)/nnod;
    
    fprintf(fid_log,' %4i | %3i x%3i x%3i | %6i (%5.1f %%) |      %2i \n',...
        itest,nsplit_now,sumSDnods(itest)-nnod,pct_shared(itest),ncomm(itest));
end

ind1      = pct_shared-min(pct_shared);
ind2      = ncomm-min(ncomm);
ind3      = 1e2*max(ind1)*invalid_split;
[~,ibest] = min(ind1+fewer_NB_wght*ind2+ind3);

nsplit       = nsplit(ibest,:);
[el2sd,topo] = assign_el2sd(GCOORD,EL2NOD,nsd,nsplit); % *SUBFUNCTION*

fprintf(fid_log,' *********************************************************\n');
fprintf(fid_log,' Will use case\n');
fprintf(fid_log,' %4i | %3i x%3i x%3i | %6i (%5.1f %%) |      %2i \n',...
            ibest,nsplit,sumSDnods(ibest)-nnod,pct_shared(ibest),ncomm(ibest));

fprintf(fid_log,'\n Subdomain topology\n');
if nsplit(1)==nsd
    fprintf(fid_log,' %2i',1:nsd);
    fprintf(fid_log,'\n xmin --> xmax\n');
elseif nsplit(2)==nsd
    fprintf(fid_log,' %2i',1:nsd);
    fprintf(fid_log,'\n ymin --> ymax\n');
elseif nsplit(3)==nsd
    fprintf(fid_log,' %2i',1:nsd);
    fprintf(fid_log,'\n zmin --> zmax\n');
else
    for iz=1:nsplit(3)
        fprintf(fid_log,'\n SDs at z-level %2i:\n',iz);
        for iy=nsplit(2):-1:1
            nn    = length( topo(:,iy,iz) );
            formt = ['  ' repmat('%2i  ',[1 nn])];
            if iy==nsplit(2) && nsplit(2)>1
                fprintf(fid_log,[formt  '  (back)\n'],topo(:,iy,iz)');
            elseif iy==1 && nsplit(2)>1
                fprintf(fid_log,[formt  '  (front)\n'],topo(:,iy,iz)');
            else
                fprintf(fid_log,[formt  '  \n'],topo(:,iy,iz)');
            end
        end
        fprintf(fid_log,' xmin --> xmax\n');
    end
end

end % END OF SUBFUNCTION find_best_split

% #########################################################################

function [el2sd,topo] = find_best_split_p(nsplit,GCOORD,EL2NOD,COMM,fewer_NB_wght,fid_log)

ntest = size(nsplit,1);
nnod  = double(max(max(EL2NOD)));
nsd   = COMM.nsd;

fprintf(fid_log,' Testing %1i configurations to find most efficient domain\n decomposition for %1i subdomains.\n',ntest,nsd);
fprintf(fid_log,' *********************************************************\n');

my_test       = COMM.myid:COMM.nsd:ntest;
n_my_test     = length(my_test);

pct_shared    = zeros(1,n_my_test);
ncomm         = zeros(1,n_my_test);
sumSDnods     = zeros(1,n_my_test);
invalid_split = zeros(1,n_my_test);

for itest=1:n_my_test
    nsplit_now = nsplit(my_test(itest),:);
    el2sd      = assign_el2sd(GCOORD,EL2NOD,nsd,nsplit_now); % *SUBFUNCTION*
    
    if length(unique(el2sd))<nsd
        invalid_split(itest)=1;
    end
    
    allSDNB = zeros(nsd-1,nsd);
    for isd=1:nsd
        NBs = el2sd( sum(ismember(EL2NOD(:,1:4),EL2NOD(el2sd==isd,1:4)),2)>0 );
        NBs = unique(NBs(NBs~=isd));
        nNB = length(NBs);
        allSDNB(1:nNB,isd) = NBs(:);
    end
    
    comm_scheme  = pairwise_comm_scheme(allSDNB);
    ncomm(itest) = size(comm_scheme,1);
    
    nnodSD = zeros(nsd,1);
    for isd=1:nsd
        nnodSD(isd) = length( unique( EL2NOD(el2sd==isd,: ) ) );
    end
    sumSDnods(itest)  = sum(nnodSD);
    pct_shared(itest) = 100*(sumSDnods(itest)-nnod)/nnod;
end

labBarrier
if COMM.myid==1
    tmp                 = pct_shared;
    pct_shared          = zeros(1,ntest);
    pct_shared(my_test) = tmp;
    tmp            = ncomm;
    ncomm          = zeros(1,ntest);
    ncomm(my_test) = tmp;
    tmp                = sumSDnods;
    sumSDnods          = zeros(1,ntest);
    sumSDnods(my_test) = tmp;
    tmp                    = invalid_split;
    invalid_split          = zeros(1,ntest);
    invalid_split(my_test) = tmp;
end
for isd=2:min(COMM.nsd,ntest)
    labBarrier
    if COMM.myid==1
        itest                = COMM.recv_1NB(isd,1);
        pct_shared(itest)    = COMM.recv_1NB(isd,2);
        ncomm(itest)         = COMM.recv_1NB(isd,3);
        sumSDnods(itest)     = COMM.recv_1NB(isd,4);
        invalid_split(itest) = COMM.recv_1NB(isd,5);
    elseif COMM.myid==isd
        COMM.send_1NB(1,my_test,1)
        COMM.send_1NB(1,pct_shared,2)
        COMM.send_1NB(1,ncomm,3)
        COMM.send_1NB(1,sumSDnods,4)
        COMM.send_1NB(1,invalid_split,5)
    end
end
labBarrier
if COMM.myid==1
    labBroadcast(1,pct_shared);
    labBroadcast(1,ncomm);
    labBroadcast(1,sumSDnods);
    labBroadcast(1,invalid_split);
else
    pct_shared    = labBroadcast(1);
    ncomm         = labBroadcast(1);
    sumSDnods     = labBroadcast(1);
    invalid_split = labBroadcast(1);
end

fprintf(fid_log,' case | # bi-sections |  # shared nodes  | # comm cycles \n');
for itest=1:ntest
    fprintf(fid_log,' %4i | %3i x%3i x%3i | %6i (%5.1f %%) |      %2i \n',...
        itest,nsplit(itest,:),sumSDnods(itest)-nnod,pct_shared(itest),ncomm(itest));
end

ind1      = pct_shared-min(pct_shared);
ind2      = ncomm-min(ncomm);
ind3      = 1e2*max(ind1)*invalid_split;
[~,ibest] = min(ind1+fewer_NB_wght*ind2+ind3);

nsplit       = nsplit(ibest,:);
[el2sd,topo] = assign_el2sd(GCOORD,EL2NOD,nsd,nsplit); % *SUBFUNCTION*

fprintf(fid_log,' *********************************************************\n');
fprintf(fid_log,' Will use case\n');
fprintf(fid_log,' %4i | %3i x%3i x%3i | %6i (%5.1f %%) |      %2i \n',...
            ibest,nsplit,sumSDnods(ibest)-nnod,pct_shared(ibest),ncomm(ibest));

fprintf(fid_log,'\n Subdomain topology\n');
if nsplit(1)==nsd
    fprintf(fid_log,' %2i',1:nsd);
    fprintf(fid_log,'\n xmin --> xmax\n');
elseif nsplit(2)==nsd
    fprintf(fid_log,' %2i',1:nsd);
    fprintf(fid_log,'\n ymin --> ymax\n');
elseif nsplit(3)==nsd
    fprintf(fid_log,' %2i',1:nsd);
    fprintf(fid_log,'\n zmin --> zmax\n');
else
    for iz=1:nsplit(3)
        fprintf(fid_log,'\n SDs at z-level %2i:\n',iz);
        for iy=nsplit(2):-1:1
            nn    = length( topo(:,iy,iz) );
            formt = ['  ' repmat('%2i  ',[1 nn])];
            if iy==nsplit(2) && nsplit(2)>1
                fprintf(fid_log,[formt  '  (back)\n'],topo(:,iy,iz)');
            elseif iy==1 && nsplit(2)>1
                fprintf(fid_log,[formt  '  (front)\n'],topo(:,iy,iz)');
            else
                fprintf(fid_log,[formt  '  \n'],topo(:,iy,iz)');
            end
        end
        fprintf(fid_log,' xmin --> xmax\n');
    end
end

end % END OF SUBFUNCTION find_best_split_p

% #########################################################################

function [el2sd,topo] = assign_el2sd(GCOORD,EL2NOD,nsd,nsplit)

elcenter = 0.25*[GCOORD(1,EL2NOD(:,1)) + GCOORD(1,EL2NOD(:,2)) + GCOORD(1,EL2NOD(:,3)) + GCOORD(1,EL2NOD(:,4));
                 GCOORD(2,EL2NOD(:,1)) + GCOORD(2,EL2NOD(:,2)) + GCOORD(2,EL2NOD(:,3)) + GCOORD(2,EL2NOD(:,4));
                 GCOORD(3,EL2NOD(:,1)) + GCOORD(3,EL2NOD(:,2)) + GCOORD(3,EL2NOD(:,3)) + GCOORD(3,EL2NOD(:,4))];

nel      = size(EL2NOD,1);
el2sd    = zeros(nel,1,'int32');

if any(nsplit==nsd)
    elSD                  = cell(nsplit);
    dir                   = find(nsplit==nsd);
    elSD_tmp              = split_mesh(elcenter(dir,:),nsplit(dir));
    for ii=1:nsd
        elSD{ii} = elSD_tmp{ii};
    end
else
    elSD = cell(nsplit);
    nx   = nsplit(1);   % number of splits in 1st dimension
    if nx>1
        dir      = 1; % direction (x-, y- or z)
        elSD_tmp = split_mesh(elcenter(dir,:),nsplit(1));
        for ix=1:nx
            elSD{ix} = elSD_tmp{ix};
        end
    else
        elSD{1} = 1:nel;
    end
    
    ny  = nsplit(2);     % number of splits in 2nd dimension
    ind = 1:nel;
    if ny>1
        dir = 2; % direction (x-, y- or z)
        for ix=1:nx
            if ~isempty(elSD{1}); ind = elSD{ix,1}; end
            elSD_tmp = split_mesh(elcenter(dir,ind),nsplit(2));
            for iy=1:ny
                els         = elSD_tmp{iy};
                elSD{ix,iy} = ind(els);
            end
        end
    end

    nz  = nsplit(3);     % number of splits in 3rd dimension
    ind = 1:nel;
    if nz>1
        dir = 3; % direction (x-, y- or z)
        for ix=1:nx
            for iy=1:ny
                if ~isempty(elSD{1}); ind = elSD{ix,iy}; end
                elSD_tmp = split_mesh(elcenter(dir,ind),nsplit(3));
                for iz=1:nz
                    els            = elSD_tmp{iz};
                    elSD{ix,iy,iz} = ind(els);
                end
            end
        end
    end
end
topo = zeros(nsplit);
for isd=1:nsd
    els        = elSD{isd};
    el2sd(els) = isd;
    topo(isd)  = isd;
end

if any(el2sd==0)
    error('Subdomain splitting failed. Look in function "assign_el2sd".')
end
if unique(el2sd)~=nsd
    error('One or more subdomains got no elements. MESH too small?')
end

end % END OF SUBFUNCTION assign_el2sd

% #########################################################################

function elSD = split_mesh(elcenter,nblocks)
    [~,iel] = sort(elcenter);
    nel     = length(iel);
    nelblo  = round(nel/nblocks);
    ii      = 1;
    jj      = nelblo;
    elSD    = cell(1,nblocks);
    for iblock=1:nblocks
        elSD{iblock} = iel(ii:jj);
        ii = jj + 1;
        if iblock==nblocks-1
            jj = nel;
        else
            jj = jj+nelblo;
        end
    end
end % END OF SUBFUNCTION split_mesh