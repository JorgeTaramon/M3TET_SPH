function [els,gX_PT,iout] = tsearch3(GCOORD,EL2NOD,gX_PT,OPTS,els)

if isempty(gX_PT)
    els = int32([]);
    return
end
if nargin<4 || ~isfield(OPTS,'fid_log')
    OPTS.fid_log = 1;
end
if nargin<5
    els = [];
    fprintf(OPTS.fid_log,' WARNING: No element list for point search provided.\n');
end
OPTS = check_options(OPTS); % *NESTED FUNCTION*

% tsearch2 becomes very slow if it searches points that are actually ouside
% the mesh. Hence points ouside the mesh must be either removed from the
% list (i.e. "nan" is returned for these points) or must be shifted onto
% the domain boundary.
[iout,gX_PT] = find_pts_outside_domain(GCOORD,gX_PT,OPTS);

EL2NOD = EL2NOD(1:4,:);
nnod   = max(EL2NOD(:));
GCOORD = GCOORD(:,1:nnod);
if exist('tsearch2','file')
    opts.nthreads = OPTS.nthreads;
    if ~isempty(els)
        els = tsearch2(GCOORD,uint32(EL2NOD([1 2 4 3],:)),gX_PT,OPTS.WS,els(:)',opts);
    else
        els = tsearch2(GCOORD,uint32(EL2NOD([1 2 4 3],:)),gX_PT,OPTS.WS,[],opts);
    end
    method = 'tsearch2';
    
else
    els    = tsearchn(GCOORD',double(EL2NOD)',gX_PT');
    method = 'tsearchn';
end

ilost      = find(isnan(els) | els==0); 
nlost      = length(ilost);
els(ilost) = 0;

if nlost==0
    if OPTS.verbose
        fprintf(OPTS.fid_log,' All %1i points have been located by %s. \n',length(els),method);
    end

else
    if OPTS.verbose
        fprintf(OPTS.fid_log,' Number of points not located by %s: %1i\n',method,nlost);
    end
    gX_PT_lost = gX_PT(:,ilost); % reduce vector sizes

    if OPTS.verbose
        fprintf(OPTS.fid_log,' Searching %1i points in nearest elements...\n',nlost);
    end
    els3       = point_search_in_nearest_tets(GCOORD,EL2NOD,gX_PT_lost,OPTS.verbose);
    els(ilost) = els3;
    
%     ilost   = find(isnan(els) | els==0); 
%     nlost   = length(ilost);
%     if nlost
%         els_bad   = els(ilost);
%         gX_PT_lost = gX_PT(:,ilost);
% 
%         ncheck = 60;   % search xx nearest elements
%         els    = find_remaining_points(els,GCOORD,EL2NOD,els_bad,gX_PT_lost,ilost,nlost,...
%             ncheck,nptblo); % *SUBFUNCTION*
%     end
end

if OPTS.check_els
    % Check that each point is in the correct element
    % ===============================================
    % Undo the permutation from line #24 (sorting points so that x-coordinate
    % increases). This is only done for the last check at the end of this fct
    eX   = local_coords_3d(GCOORD,EL2NOD,els,gX_PT);
    tol  = 1e-4;
    iok  = all(eX>=0-tol,1) & all(eX<=1+tol,1) & sum(eX)<=1+tol;
    ibad = find(~iok);
    if ~isempty(ibad)
        fprintf(OPTS.fid_log,' %1i points located in wrong elements...',length(ibad));
    else
        fprintf(OPTS.fid_log,' Checked "lc", all points located in correct elements.\n');
    end
end

% =========================================================================
%                            NESTED FUNCTIONS
% =========================================================================

function OPTS = check_options(OPTS)

% DEFAULTS (IF NOT PROVIDED)
if ~isfield(OPTS,'fid_log')
    OPTS.fid_log = 1;
end
if ~isfield(OPTS,'ztol')
    OPTS.ztol = 1e-6 * (max(GCOORD(3,:))-min(GCOORD(3,:)));
end
if ~isfield(OPTS,'verbose')
    OPTS.verbose = 0;
end
if ~isfield(OPTS,'check_els')
    OPTS.check_els = 0;
end
if ~isfield(OPTS,'WS')
    OPTS.WS.NEIGHBORS = calc_tetra_neighbors(EL2NOD);
    fprintf(OPTS.fid_log,' WARNING: Have to calculate element neighbor information. This slows down tsearch3!\n');
end
if ~isfield(OPTS,'nthreads')
    OPTS.nthreads = 0;
end
if ~isfield(OPTS,'EL2NOD_top')
    OPTS.EL2NOD_top = [];
%     fprintf(OPTS.fid_log,' WARNING: No information on top boundary provided. This slows down tsearch3!\n');
end
if ~isfield(OPTS,'EL2NOD_bot')
    OPTS.EL2NOD_bot = [];
%     fprintf(OPTS.fid_log,' WARNING: No information on bottom boundary provided. This slows down tsearch3!\n');
end
if ~isfield(OPTS,'xmin')
    OPTS.xmin = min(GCOORD(1,:));
end
if ~isfield(OPTS,'xmax')
    OPTS.xmax = max(GCOORD(1,:));
end
if ~isfield(OPTS,'ymin')
    OPTS.ymin = min(GCOORD(2,:));
end
if ~isfield(OPTS,'ymax')
    OPTS.ymax = max(GCOORD(2,:));
end
if ~isfield(OPTS,'zmin')
    OPTS.zmin = min(GCOORD(3,:));
end
if ~isfield(OPTS,'zmax')
    OPTS.zmax = max(GCOORD(3,:));
end

end % END OF NESTED FUNCTION check_options

end % END OF FUNCTION tsearch3

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function els = point_search_in_nearest_tets(GCOORD,EL2NOD,gX_PT,verbose)

if nargin==0
    tmp    = load('TestData_point_search');
    GCOORD   = tmp.GCOORD;
    EL2NOD = tmp.EL2NOD;
    gX_PT  = tmp.gX_PT_bad;
    clear tmp
end

FigNo     = 0;
nptblo    = 5000;  % number of points handled at once
tol       = 1e-2;  % m (to pre-select possible elements using xmin, xmax, ...)
ncheck    = [];    % maximum number of elements that will be tested
                   % leave empty ([]) to test all possible elements

gX_PT(1,:) = min(gX_PT(1,:),max(GCOORD(1,:)));
gX_PT(1,:) = max(gX_PT(1,:),min(GCOORD(1,:)));
gX_PT(2,:) = min(gX_PT(2,:),max(GCOORD(2,:)));
gX_PT(2,:) = max(gX_PT(2,:),min(GCOORD(2,:)));
gX_PT(3,:) = min(gX_PT(3,:),max(GCOORD(3,:)));
gX_PT(3,:) = max(gX_PT(3,:),min(GCOORD(3,:)));

[npe,nel] = size(EL2NOD);
nPT       = size(gX_PT,2);
els       = nan(nPT,1);
nptblo    = min(nPT,nptblo);
nblo      = ceil(nPT/nptblo);

% sort lost points so that their x-coordinate increases
[~,perm_PT] = sort(gX_PT(1,:));
gX_PT       = gX_PT(:,perm_PT);
% perm_PT = 1:nPT;
% 
% points are like element centers that need to be sorted into subdomains
% except that the number of subdomains is undefined and you have to create SDs
% until the number of points in each SD is below a threshold
% see "craete_subdomains"

elc       = calc_tetra_center(GCOORD,EL2NOD);
x_el      = reshape(GCOORD(1,EL2NOD(:,1:nel)),npe,nel);
y_el      = reshape(GCOORD(2,EL2NOD(:,1:nel)),npe,nel);
z_el      = reshape(GCOORD(3,EL2NOD(:,1:nel)),npe,nel);

il         = 1;
els1       = zeros(1,nPT*40,'int32'); % estimating a maximum of 40...
gX_PT2     = zeros(3,nPT*40);         % ...possible elements per point
ptr_2_PT   = zeros(1,nPT*40,'int32'); %
ipt        = 1;
jpt        = nptblo;
nmax       = 0;
for iblo = 1:nblo
    gX_PT_BLO   = gX_PT(:,ipt:jpt);
    perm_PT_blo = perm_PT(ipt:jpt);
    els_imposs  = find( min(x_el)-tol > max(gX_PT_BLO(1,:)) | ...
                        max(x_el)+tol < min(gX_PT_BLO(1,:)) | ...
                        min(y_el)-tol > max(gX_PT_BLO(2,:)) | ...
                        max(y_el)+tol < min(gX_PT_BLO(2,:)) | ...
                        min(z_el)-tol > max(gX_PT_BLO(3,:)) | ...
                        max(z_el)+tol < min(gX_PT_BLO(3,:)) );
    ptr_2_allels = setdiff(1:nel,els_imposs);
    
    xmin_el = min(x_el(:,ptr_2_allels))-tol;
    xmax_el = max(x_el(:,ptr_2_allels))+tol;
    ymin_el = min(y_el(:,ptr_2_allels))-tol;
    ymax_el = max(y_el(:,ptr_2_allels))+tol;
    zmin_el = min(z_el(:,ptr_2_allels))-tol;
    zmax_el = max(z_el(:,ptr_2_allels))+tol;

    if FigNo
        if iblo==1
            cols = jet(nblo);
            figure(FigNo);clf;set(gcf,'Renderer','OpenGL');
            plot(xz_top(1,:),xz_top(2,:),'r-'); hold on
            plot(xz_bot(1,:),xz_bot(2,:),'r-');
        end
        scatter(gX_PT_BLO(1,:)',gX_PT_BLO(3,:)',10,cols(iblo,:));
    end

    for i=1:nptblo % loop over lost points
        % find elements that can possibly contain the point:
        els_possib = find(       xmin_el             <= gX_PT_BLO(1,i) );
        els_possib = els_possib( xmax_el(els_possib) >= gX_PT_BLO(1,i) );
        els_possib = els_possib( ymin_el(els_possib) <= gX_PT_BLO(2,i) );
        els_possib = els_possib( ymax_el(els_possib) >= gX_PT_BLO(2,i) );
        els_possib = els_possib( zmin_el(els_possib) <= gX_PT_BLO(3,i) );
        els_possib = els_possib( zmax_el(els_possib) >= gX_PT_BLO(3,i) );
        els_possib = ptr_2_allels( els_possib );
        if any(els_possib==0)
            error('This should not happen.');
        end
        n = length(els_possib);
        if ~n
            error('This should not happen.');
        end

        if ~isempty(ncheck)
            % sort the elements so that distance between element center and the
            % lost point increases
            dist       = sum( (repmat(gX_PT_BLO(:,i),1,n) - elc(:,els_possib)).^2 ); % sqrt skipped for speed (does not affect sort(dist))
            [~,p]      = sort(dist);
            els_possib = els_possib(p);
            n          = min(ncheck,length(els_possib));
        else
            n          = length(els_possib);
        end
        nmax              = max(nmax,n);
        iu                = il+n-1;
        els1(il:iu)       = els_possib(1:n);
        gX_PT2(:,il:iu)   = gX_PT_BLO(:,i) * ones(1,n);
        ptr_2_PT(il:iu)   = int32(perm_PT_blo(i));
        il                = iu + 1;
    end
    ipt = jpt+1;
    jpt = ipt+nptblo-1;
    if iblo==nblo-1
        jpt    = min(jpt,nPT);
        nptblo = jpt-ipt+1;
    end
end
clear gX_PT_BLO els_possib els_imposs
els1(il:end)     = [];
gX_PT2(:,il:end) = [];
ptr_2_PT(il:end) = [];

% Calculate local coordinates for all points in all possible elements
tols     = [0 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1];
eX       = local_coords_3d(GCOORD,EL2NOD,els1,gX_PT2);
iPT_noel = 1:nPT;
iok      = false(size(ptr_2_PT));
for it=1:length(tols)
    tol = tols(it);
    % check which local coords a in a valid range
    ptr_noel = find(ismember(ptr_2_PT,iPT_noel));
    iok_itol = false(size(ptr_2_PT));
    
    iok_itol(ptr_noel) = all(eX(:,ptr_noel)>=0-tol,1) & ...
                         all(eX(:,ptr_noel)<=1+tol,1) & ...
                         sum(eX(:,ptr_noel))<=1+tol;
           
    iok      = iok | iok_itol;
    
    iPT_noel = setdiff(1:nPT,unique(ptr_2_PT(iok)));
    nok      = length(iPT_noel);
    if nok
        if verbose
            fprintf(' LCOORD-tolerance %.0e: %1i points have no element.\n',tol,nok);
        end
    else
        if verbose
            fprintf(' All points have at least one possible element.\n');
        end
        break
    end
end

% Reduce lists to the possible elements
els1        = els1(iok);
ptr_2_PT    = ptr_2_PT(iok);
eX          = eX(:,iok);

% Find the points for which only ONE element is possible
[~,ia,ib]    = unique(ptr_2_PT);
nel_possib   = accumarray(ib(:),ones(size(ib(:)))); clear ib
iok_1el      = find(nel_possib==1);
els(iok_1el) = els1(ia(iok_1el));

if any(nel_possib>1)
    ptr_2_PT(ia(iok_1el)) = [];
    els1(ia(iok_1el))     = [];
    eX(:,ia(iok_1el))     = [];

    % Pick the element with the "better" local coords.
    rstu     = [eX; 1-sum(eX,1)]; clear eX
    sum_err  = sum( abs(min(0,rstu)) ) + sum( abs(min(0,1-rstu)) );
    [~,perm] = sort(sum_err,'descend'); % descend to put "best" elements in last position
    
    % Because the elements are sorted with the "best" element per point in
    % the last position, the following will write the "best" element for
    % each point into els.
    ptr_2_PT      = ptr_2_PT(perm);
    % a([1 1 2 2 2]) = 1:5;
    %   produces a = [2 5]
    els(ptr_2_PT) = els1(perm);
end
nlost = length(find(isnan(els)));

if ~nlost
    if verbose
        fprintf(' All %1i points have been located.\n',nPT);
    end
else
    fprintf(' Number of not located points: %1i\n',nlost);
    error(' This must not happen. Check tolerances etc.');
    if FigNo
        gX_PT(:,igood) = [];
        scatter(gX_PT(1,ibad)',gX_PT(3,ibad)',20,'k');
    end
%         filename = ['Lab' num2str_d(numlabs,2) 'x' num2str_d(labindex,2) '_' num2str(nPT) 'lostPTs.dat'];
%         fid      = fopen(filename,'w');
%         fprintf(fid,'%10.3f %10.3f %10.3f\n',gX_PT);
%         fclose(fid);
end

end % END OF SUBFUNCTION point_search_in_nearest_tets

% #########################################################################

% function ptr_2_grp = group_all_points(gX_PT,max_grp_size)
% 
% nPT  = size(gX_PT,2);
% ngrp = ceil(nPT/max_grp_size);
% while 1
% 	n = log2(ngrp);
%     if n==round(n)
%         break
%     end
%     ngrp = ngrp + 1;
% end
% 
% iexp   = 1:log2(ngrp);
% k      = [1 2.^iexp];
% nsplit = [k(end:-1:1)' k' ones(ngrp+1,1)];
% 
% end % END OF SUBFUNCTION group_all_points

% #########################################################################
% #########################################################################

% function el1 = check_lcoords(GCOORD,EL2NOD,els1,gX_PT)
%     tol  = 1e-8;
%     if size(gX_PT,2)==1
%         gX_PT = repmat(gX_PT,1,length(els1));
%     end
%     eX   = local_coords_3d(GCOORD,EL2NOD,els1,gX_PT);
%     iok  = find(all(eX>=0-tol,1) & all(eX<=1+tol,1) & sum(eX)<=1+tol);
%     if ~isempty(iok)
%         el1 = els1(iok(1));
%     else
%         el1 = 0;
%     end
% end

% #########################################################################

% other algorithms...

% case 0
%     nel       = size(EL2NOD,1);
%     nel_check = min(nel,10);
%     if verbose
%         fprintf(' Fix 2: Check local coords in %1i els close to each point.\n',nel_check);
%     end
%     elc    = calc_tetra_center(GCOORD',EL2NOD');
%     for i=1:nPT
%         dist       = sum( (repmat(gX_PT_bad(i,:)',1,nel) - elc).^2 ); % sqrt skipped for speed (does not affect sort(dist))
%         [~,els1]   = sort(dist);
%         els_bad(i) = check_lcoords(GCOORD',EL2NOD',els1(1:nel_check),gX_PT_bad(i,:)'); % *NSUBFUNCTION*
%     end
%     
% case 1
%     if verbose
%         fprintf(' Fix 1: Calling tsearchn with scaled coordinates.\n');
%     end
%     Lxyz    = max(GCOORD)-min(GCOORD);
%     scale   = 1./mean(Lxyz);
%     els_bad = tsearchn(scale*GCOORD,double(EL2NOD),scale*gX_PT_bad);

% case 2
%     if verbose
%         fprintf(' Fix 3: Check local coords in els connected to nearest node.\n');
%     end
%     nearest_node = dsearchn(GCOORD,double(EL2NOD),gX_PT_bad);
%     for i=1:nPT
%         els1         = find( sum(ismember(EL2NOD,nearest_node(i)),2) );
%         els(ibad(i)) = check_lcoords(GCOORD',EL2NOD',els1,gX_PT_bad(i,:)'); % *NSUBFUNCTION*
%     end
% 
% case 3
%     if verbose
%         fprintf(' Fix 3: Slightly change point coords and search again.\n');
%     end
%     nel   = size(EL2NOD,1);
%     shake = [1e-4 1e-4 1e-4 1e-4...
%              1e-3 1e-3 1e-3 1e-3...
%              1e-2 1e-2 1e-2 1e-2] * min(Lxyz)./nel^(1/3);
%     for i=1:length(shake)
%         gX_PT_bad  = gX_PT(ibad,:) + shake(i)*rand(nPT,3);
%         els(ibad)  = tsearchn(GCOORD,double(EL2NOD),gX_PT_bad);
%         ibad       = isnan(els) | els==0; 
%         nPT       = length(find(ibad));
%         if ~nPT
%             break
%         end
%     end
% 
% case 4
%     if verbose
%         fprintf(' Fix 4: Check local coords in all els of the mesh.\n');
%     end
%     nel    = size(EL2NOD,1);
%     nelblo = 50000;
%     il     = 1;
%     iu     = min(nel,nelblo);
%     nblo   = ceil(nel/nelblo);
%     for iblo=1:nblo
%         for i=1:nPT
%             eX  = local_coords_3d(GCOORD',EL2NOD',(il:iu),gX_PT(ibad(i),:)');
%             iok = find(all(eX>=0,1) & all(eX>=1,1) & sum(eX)<=1);
%             if ~isempty(iok)
%                 els(ibad(i)) = els1(iok);
%             end
%         end
%         ibad = find(isnan(els) | els==0); 
%         nPT = length(ibad);
%     end