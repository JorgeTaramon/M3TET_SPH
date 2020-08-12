function [vars_PT,els,locinfo] = locate_points_interp_vars_3d_p_sph...
    (MESH,COMM,gTH_PT,vars,OPTS_search,OPTS_interp,els)
% Usage: [vars_PT,els,locinfo] = locate_points_interp_vars_3d_p_sph...
%   (MESH,COMM,gX_PT,vars,OPTS_search,OPTS_interp,els)
%
% Purpose: Locate points with coordinates gX_PT in the mesh and
%          interpolates variables at these locations
%
% Input:
%   MESH     : [structure]    : FE mesh parameters
%   COMM     : [structure]    : inter-subdomain communication data
%   gX_PT    : [matrix]       : coordinates of points to be located (3 x nnod)
%   vars     : [matrix]       : variables to be interpolated (nnod x nvar)
%   OPTS_search : [structure] : options for tsearch2 (e.g.
%   OPTS_interp : [structure] : options for interpolation (e.g. method)
%   els      : [vector]       : guess for which element contains each point
% Output:
%   vars_PT  : [matrix]    : variables interpolated at points (npt x nvar)
%   els      : [vector]    : elements in which each point was located
%   locinfo  : [structure] : profiling data
%
% written by J.Hasenclever, 2011
% Email contact: jhasenclever@geomar.de
%
% For numerical methods see online Ph.D. thesis
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)

% JH Sep 2012
% JH Dec 2014 : added tsearch2 and tested for fastest method
% JMT Nov 2016 : Now it works for spherical coordinates although only for
%                COMM.nsd = 1 and nmg = 1 

warning('locate_points_interp_vars_3d_p_sph is only working for COMM.nsd = 1 and nmg = 1')

GCOORD     = MESH.GCOORD;
GCOORD_SPH = MESH.GCOORD_SPH;
EL2NOD     = MESH.EL2NOD;
nmg        = length(MESH.EL2NOD);
if COMM.nsd>1
    GCOORD_D = MESH.GCOORD_D;
    EL2NOD_D = MESH.EL2NOD_D;
end
nnod       = size(GCOORD_SPH,2);
if size(vars,2)==nnod; vars=vars'; end
npt        = size(gTH_PT,2);
nvar       = size(vars,2);
tol_lc     = 5e-4; % by default it was 1e-4 but for some nodes located  
                   % inside isoparametric elements (curved elements) this 
                   % tol_lc was not enough. For tol_lc = 5e-3 it works fine
                   % FIXME
% Profiling data
locinfo.npt    = npt;
locinfo.npt_lc = zeros(1,nmg);
locinfo.t_lc   = zeros(1,nmg);

% If no (useful) vector "els" is provided: Calculate a guess for the
% element in which each each point is located
% ==================================================================
if nargin<7 || length(els)~=npt
    if npt==nnod
        % Assuming PTs are backtracking points
        % --> choose one element connected to each node
        els      = guess_elems(EL2NOD{1}); % *SUBFUNCTION*
        check_lc = 0; % 1 --> calc local coords; find points to be located
    else
        % No idea how to construct a good guess..
        check_lc = 0; % 1 --> calc local coords; find points to be located
    end

else 
    check_lc = 0; % 1 --> calc local coords; find points to be located
    ibad = find(els==0 | isnan(els));
    if ~isempty(ibad)
        if npt==nnod
            el_guess  = guess_elems(EL2NOD{1}); % *SUBFUNCTION*
            els(ibad) = el_guess(ibad);
            if length(ibad)/npt > 0.1
                check_lc = 0; % 1 --> calc local coords; find points to be located
            end
        else
            els(ibad) = uint32(1);
            check_lc  = 0; % 1 --> calc local coords; find points to be located
        end
    end
end

img_loc = zeros(1,npt);
if check_lc
    % Check local coordinates assuming points are in els
    lc          = local_coords_3d_sph(GCOORD_SPH,EL2NOD{1},els,gTH_PT,MESH);
    isrch       = any(lc>1,1) | any(lc<0,1) | sum(lc,1)>1;
    
    % Points with index "isrch" will be searched on the 
    % coarsest mesh and then recursively upwards
    % All others have already correct "els" and "lc"
    lc(:,isrch) = 0;
    
    % This marks the points already located (see line beginning with
    % "ipt_in = find...")
    img_loc(~isrch) = -1;
    isrch           = find(isrch);
    
    % Calculate grand-parent elements on coarsest mesh for the points that are
    % in the wrong element  on the fine mesh
    els(isrch) = uint32(ceil( double(els(isrch))./(8^(nmg-1)) ));
else
    isrch = 1:npt;
    lc    = zeros(3,npt);
    els   = ones(npt,1,'uint32');
end

tic
if COMM.nsd>1
    error('it needs to be coded')
    % Search points in coarse domain mesh
    [els_D,gX_PT(:,isrch)] = tsearch3...
        (GCOORD_D,EL2NOD_D(1:4,:),gX_PT(:,isrch),OPTS_search,els(isrch));
    locinfo.t_srch         = toc;
    locinfo.npt_srch       = length(isrch);tic
    iloc                   = els_D>0;
    % check in which subdomain each element of list "els_D" lies
    iSDel                  = zeros(size(els_D));
    iSDel(iloc)            = MESH.el2sd(els_D(iloc));
    % find those that are within SD==myid
    myels                  = iSDel==COMM.myid;
    % convert element numbers from domain to subdomain numbering;
    % write these "subdomain" elements into the list of elements
    els(isrch(myels))      = MESH.el_D2SDc( els_D(myels) );
    % update the multigrid level list
    img_loc(isrch(myels))  = nmg; % points that were located within the SD
    
    % PARALLEL SECTION: Some points have been located n oterh subdomains
    % and this SD has located points from other SDs. These points must be
    % communicated so that every SD has a complete list of points.
    % =====================================================================
    [els_D_NB,gX_PT_NB,ptr2nbPTs,ptr2myPTs,NBs] = communicate_pts...
        (COMM,els_D,iSDel,gX_PT,isrch,OPTS_search.verbose); % *SUBFUNCTION*

    % Merge the list of points inside this SD and the list of points that have
    % to be located for neighboring SDs
    % =========================================================================
    nSDel_all = length(els_D_NB);
    if nSDel_all
        img_loc = [img_loc nmg*ones(1,nSDel_all)];
        els     = [els(:); MESH.el_D2SDc(els_D_NB)']; % GLOBAL element # --> SD element #
        gX_PT   = [gX_PT gX_PT_NB];
    end
    locinfo.t_comm = toc;
else
    % Search points in coarse mesh
%     [els(isrch),gX_PT(:,isrch)] = tsearch3...
%         (GCOORD,EL2NOD{nmg}(1:4,:),gX_PT(:,isrch),OPTS_search,els(isrch));

%     [els(isrch),gTH_PT(:,isrch),OPTS_search] = tsearch2_sph_old_version...
%         (GCOORD_SPH,EL2NOD{1},gTH_PT(:,isrch),OPTS_search,MESH);

    [els(isrch),gTH_PT(:,isrch),OPTS_search] = tsearch2_sph...
        (GCOORD_SPH,EL2NOD{1},gTH_PT(:,isrch),OPTS_search,MESH);
    locinfo.t_srch       = toc;
    locinfo.npt_srch     = length(isrch);
    iloc                 = els(isrch)>0;
    img_loc(isrch(iloc)) = nmg; % points that were located within the SD
end


% Loop from coarse to the finest geometric MG mesh and calculate the child 
% element in which the points are located. I.e.: based on the local
% coordinates of a point in an element it is calculated in which of the 8
% subelements the point is located on the next finer MG level. This info is
% used to calculate the element number on the next finer mesh.
% =======================================================================
for img=nmg:-1:1
                                                                            tic
    ipt       = find(img_loc==img);
    lc(:,ipt) = local_coords_3d_sph(GCOORD_SPH,EL2NOD{1},els(ipt),gTH_PT(:,ipt),MESH);
    inot_loc  = any(lc(:,ipt)<0-tol_lc) | any(lc(:,ipt)>1+tol_lc) | ...
                sum(lc(:,ipt))>1+tol_lc;
    inot_loc  = ipt(inot_loc);
    if ~isempty(inot_loc)
        inot_loc      = ipt(inot_loc);
        [els(inot_loc),lc(:,inot_loc)] = locate_points_in_neighboring_element...
            (els(inot_loc),gTH_PT(:,inot_loc),GCOORD_SPH,EL2NOD{1},MESH);
%         % Get  point with worst local coordinates
%         [~,jpt] = max( max(abs(lc(:,inot_loc)))  );
%         jpt     = inot_loc(jpt);
%         lc(:,jpt)
%         el      = els(jpt);
%         % Find all elements connected to element "jel"
%         el_nb   = find(any(ismember(EL2NOD{1},EL2NOD{1}(1:4,el)),1));
%         nel_nb  = length(el_nb);
%         % Assume the point is in the element neighbors and calculate local
%         % coordinates
%         lc_nb   = local_coords_3d_sph(GCOORD_SPH,EL2NOD{1},el_nb,...
%             repmat(gTH_PT(:,jpt),1,nel_nb),MESH);
%         lc_nb(4,:) = 1-sum(lc_nb(1:3,:));
%         fprintf('el=%6i lc=%+.2e %+.2e %+.2e\n',el,lc(:,jpt));
%         fprintf('el=%6i lc=%+.2e %+.2e %+.2e %+.2e\n',[el_nb(:) lc_nb']');
%         el_ok = el_nb(all(lc_nb>=0 & lc_nb<=1,1));
%         if isempty(el_ok)
%             lc_all      = local_coords_3d_sph(GCOORD_SPH,EL2NOD{1},1:MESH.nel,...
%                 repmat(gTH_PT(:,jpt),1,MESH.nel),MESH);
%             lc_all(4,:) = 1-sum(lc_all(1:3,:));
%             el_ok       = find(all(lc_all>=0 & lc_all<=1,1));
%         end
%         % save(['Workspace_Lab' COMM.prefix '_location_error'])
%         error(' This should not happen: %1i points on MG%1i wrong local coords!!',...
%             length(find(inot_loc)),img);
    end
    
    % Adjust local coordinates so that  0 <= [r,s,t,u=1-r-s-t] <= 1
    lc(:,ipt)    = adjust_local_coords( lc(:,ipt) ); % *SUBFUNCTION*
    
%     % Recalculate the physical coordinate using the adjusted local
%     % coordinates to avoid problems in next iteration.
%     % Use LINEAR interpolation here (only use vertex nodes EL2NOD(1:4,:)!
%     xyz_interp   = interp3d_tet(EL2NOD{img}(1:4,:),els(ipt),lc(:,ipt),GCOORD_SD');
%     gX_PT(:,ipt) = xyz_interp';
    
    if img==1                                                       
        break
    end
    
    % As long as we're not on the finest level:
    % Use els and lc to calculate the child element within
    % which gX_PT is located 
    els(ipt)     = calc_tetra_childs(els(ipt),lc(:,ipt));
    img_loc(ipt) = img-1;
                                                                            locinfo.t_lc(img)=locinfo.t_lc(img)+toc;
end
locinfo.t_lc(img)=locinfo.t_lc(img)+toc;
clear EL2NOD

% if ~check_lc
%     save('els_no_lc_check','els','lc');
% else
%     data = load('els_no_lc_check');
%     max(abs(els-data.els))
%     max(abs(lc-data.lc)')
% end

if isfield(OPTS_search,'verbose') && OPTS_search.verbose
    fprintf(' Number of points found in old element   : %1i (%5.2f%%)\n',...
            npt-locinfo.npt_srch,100*(npt-locinfo.npt_srch)/npt);
    fprintf(' Number of points searched in coarse mesh: %1i (%5.2f%%)\n',...
            locinfo.npt_srch,100*locinfo.npt_srch/npt);
    fprintf(' Time for searching on coarse mesh       : %5.2f\n',locinfo.t_srch);
    fprintf(' Time for re-locating on finer meshes    :');
    fprintf(' %5.2f ',locinfo.t_lc);
    fprintf(' \n');
end

% if numlabs>1
%     filename = ['Workspace_Lab' num2str_d(labindex,2)]
%     save(filename);
%     stop
% else
%     myid     = 1;
%     filename = ['Workspace_Lab' num2str_d(myid,2)]
%     load(filename);
% end
    
% All points are located on the finest MG level. Now interpolate the
% variable(s) at points that this SD needs (ipt_in) and that other SDs need
% (ipt_4NB) but lie within this SD.
% =========================================================================
ipt_in  = find( img_loc(1:npt)); % index of points located inside this SD
ipt_4NB = npt+1:length(img_loc); % positions where points for other SDs are
                                 % stored
vars_PT = zeros(npt,nvar);
npt_in  = length(ipt_in);

% Interpolate variables at points
% ===============================
ind_els = [ipt_in ipt_4NB];
nnodel  = size(MESH.EL2NOD{1},1);
if nnodel==10 && ~strcmp(OPTS_interp.method_interp,'quadratic')
    EL2NOD = tetmesh_p2_to_p1(GCOORD,MESH.EL2NOD{1});
    error('it needs to be debugged')
    %==========================================================================================================
    % Compute elements in relation with the cone and crossing phi = 2pi for the new connectivity 
    % For a 4-node connectivity there are no isoparametric elements (curved-edge elements). 
    % Elements are classified in:
    %   - els_out_cone_no_cross_2pi
    %   - els_out_cone_cross_2pi
    %   - els_in_cone_no_iso
    MESH.els_in_cone_iso           = [];
    [MESH.els_in_cone_no_iso,~,~,els_out_cone,~] = check_els_in_cone(GCOORD_SPH,EL2NOD,MESH.theta_cone);
    [els_cross_2pi,~]              = check_phi(GCOORD_SPH,EL2NOD(:,els_out_cone));
    MESH.els_out_cone_cross_2pi    = els_out_cone(els_cross_2pi);
    MESH.els_out_cone_no_cross_2pi = els_out_cone;
    MESH.els_out_cone_no_cross_2pi(ismember(MESH.els_out_cone_no_cross_2pi,MESH.els_out_cone_cross_2pi)) = [];
    %==========================================================================================================
    els2   = calc_tetra_childs(els(ind_els),lc(:,ind_els));
    lc2    = local_coords_3d_sph(GCOORD_SPH,EL2NOD,els2,gTH_PT(:,ind_els),MESH);
else
    EL2NOD = MESH.EL2NOD{1}; % 10-node connectivty matrix
    els2   = els(ind_els);
    lc2    = lc(:,ind_els);
end
if strcmp(OPTS_interp.method_interp,'cubic')
    tmp = interp3d_cubic_p(GCOORD,EL2NOD,COMM,els2,lc2,vars,OPTS_interp); % interpolate
else
    tmp = interp3d_tet(EL2NOD,els2,lc2,vars,OPTS_interp); % interpolate
end
els(npt+1:end) = [];

if npt_in
    vars_PT(ipt_in,:) = tmp(1:npt_in,:);
end

if COMM.nsd>1
    npt_4NB = length(ipt_4NB);
    if npt_4NB
        vars_PT_for_NB = tmp(npt_in+1:end,:);
    else
        vars_PT_for_NB = [];
    end
    clear tmp

    if OPTS_interp.verbose && npt_4NB
        fprintf('\n Interpolation included %1i pts for other SDs. Now sending results.\n',npt_4NB);
    end
    % Communicate velocities that are needed by other subdomains (i.e. BTs that
    % are located in this SD but belong to nodes of other SDs)
    % =========================================================================
    for iNB=1:length(NBs)
        NB = NBs(iNB);
        labBarrier
        if NB>0
            iloc  = ptr2nbPTs==NB;                                          locinfo.npt_for_SD(NB)=length(find(iloc));
            if OPTS_interp.verbose
                fprintf(' SD%1i: sending %1i and receiving %1i values\n',...
                         NB,length(find(iloc)),length(ptr2myPTs{NB})); end
            sdata = vars_PT_for_NB(iloc,:);
            rdata = COMM.sendnrecv_1NB(NB,sdata,103);                       locinfo.npt_from_SD(NB)=length(ptr2myPTs{NB});
            if length(ptr2myPTs{NB})~=size(rdata,1);
                error(' Received %1i values instead of the expected %1i from SD%1i',...
                      length(rdata),length(ptr2myPTs{NB}),NB);
%             else
%                 if OPTS_interp.verbose; fprintf(' Received the correct number of velocities\n'); end
            end
            if ~isempty(rdata)
                vars_PT(ptr2myPTs{NB},:) = rdata;
%                 els(ptr2myPTs{NB})       = el_guess(ptr2myPTs{NB});
            end
%         else
%             if OPTS_interp.verbose; fprintf('  Have to wait during this comm cycle...\n'); end
        end
    end
end

end % END OF FUNCTION locate_points_interp_vars_3d_p

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function lc = adjust_local_coords(lc)
% Make sure all local cooords are in the valid range 0>=lc>=1
lc( lc<0 ) = 0;
lc( lc>1 ) = 1;
ibad       = sum(lc)>1;
if any(ibad)
    scale      = repmat(sum(lc(:,ibad)),3,1);
    lc(:,ibad) = lc(:,ibad)./scale;
end
% % Make sure all local cooords are in the valid range 0>=lc>=1
% % In addition: slightly shift all points towards the element center to
% %              avoid round-off errors when claculating sub-elements
% tol        = 1e-8; % amount by which all points are shifted towards element centers
% lc         = max(tol,lc);
% lc         = min(1-tol,lc);
% scale      = repmat(sum(lc),3,1);
% scale      = scale + tol;
% lc         = lc./scale;
end % END OF SUBFUNCTION adjust_local_coords

% #########################################################################
  
function [els_D_NB,gX_PT_NB,ptr2nbPTs,ptr2myPTs,NBs] = communicate_pts...
             (COMM,els_D,iSDel,gX_PT,isrch,verbose)

myid = COMM.myid;
nsd  = COMM.nsd;

% Send points outside of this SD to the SD that includes the point.
% Receive points from other SDs that are located inside this SD.

% There are 2 ways of doing this communication:
% 1) Every SD talks to all other SDs
%    E.g.: communication for a 4-SD configuration
%    1-->2, 2-->1
%    1-->3, 3-->1
%    1-->4, 4-->1
%    2-->3, 3-->2
%    2-->4, 4-->2
%    3-->4, 4-->3
% 2) Find out which SDs have to talk and create a communication scheme that
%    allows pairwise communication without deadlock.
% Version 2 is better as it reduces the amount of comminication, especially
% for a larger number of subdomains.
nPT_from_SDi = zeros(1,nsd);
nPT_to_SDi   = zeros(1,nsd);
for isd=1:nsd
    if isd==myid
        continue
    end
    % number of PTs that will be sent to SD "isd"
    nPT_to_SDi(isd)   = length(find(iSDel==isd));
    
    % nPT_from_SDi(isd) is the number of points that SD "isd" will send
    nPT_from_SDi(isd) = COMM.sendnrecv_1NB(isd,nPT_to_SDi(isd));
    % if verbose; fprintf(' Will send %1i pts to SD%1i and will receive %1i pts\n',nPT_to_SDi(isd),isd,nPT_from_SDi(isd) ); end
end

if verbose
    fprintf(' Number of points that SD%1i will send to SD 1,2,etc :\n',myid);
    fprintf(' %1i',nPT_to_SDi);
    fprintf('\n');
    fprintf(' Number of points that SD%1i will receive from SD 1,2,etc:\n',myid);
    fprintf(' %1i',nPT_from_SDi);
    fprintf('\n');
end

% Lab 1 will now receive all required communication among SDs, then
% calculate an efficient pairwise communication scheme and broadcast.
SD2talk2 = find(nPT_to_SDi | nPT_from_SDi);
allSDNB  = zeros(nsd-1,nsd);
allSDNB(1:length(SD2talk2),1) = SD2talk2;
if myid==1
    for isd=2:nsd
        % SD "1" receives list from SD "isd"
        SD2talk2 = COMM.recv_1NB(isd); % list
        allSDNB(1:length(SD2talk2),isd) = SD2talk2(:); % save list in matrix
    end
else
    % list of SDs that SD "isd" needs to talk to
    SD2talk2 = find(nPT_to_SDi | nPT_from_SDi); 
    % SD "isd" sends list to SD 1
    COMM.send_1NB(1,SD2talk2); %#ok<FNDSB>
end
labBarrier

% Calculate pairwise, non-locking communication scheme and broadcast
if myid==1
    comm_scheme = pairwise_comm_scheme(allSDNB);
    if verbose; fprintf(' SD1 calculated this communication scheme\n');
        fprintf([repmat('%2i ',1,COMM.nsd) '\n'],comm_scheme'); 
    end
    labBroadcast( 1 , comm_scheme );
else
    comm_scheme = labBroadcast( 1 );
    if verbose; fprintf(' Received this communication scheme from SD1\n');
        fprintf([repmat(' %2i',1,COMM.nsd) '\n'],comm_scheme'); 
    end
end

% column myid is comm scheme for this SD
% rows contain the communication partner in each cycle
NBs       = comm_scheme(:,myid)';
nSDel_all = sum(nPT_from_SDi);  % total number of points this SD will receive
gX_PT_NB  = zeros(3,nSDel_all);
els_D_NB  = zeros(1,nSDel_all);
ptr2nbPTs = zeros(1,nSDel_all);
ptr2myPTs = cell(1,length(NBs));
i         = 1;
for iNB=1:length(NBs)
    NB = NBs(iNB);
    labBarrier
    if NB>0
        if verbose; fprintf(' Sending %1i pts to SD%1i, will receive %1i pts.\n',nPT_to_SDi(NB),NB,nPT_from_SDi(NB)); end
        iloc        = iSDel==NB;
        elsD_for_NB = els_D(iloc);
        ipt_for_NB  = isrch(iloc);
         % Note: elsD_for_NB is the GLOBAL element number. The next line
         % sends elsD_for_NB and receives the neighbor's elsD_for_NB list.
        rdata1 = COMM.sendnrecv_1NB(NB,elsD_for_NB,101);
        sdata  = gX_PT(:,ipt_for_NB);
        rdata2 = COMM.sendnrecv_1NB(NB,sdata,102);

        % Save this pointer to put the velocities received later from NB
        % into the correct poisiton of the velocity array of this SD
        ptr2myPTs{NB} = ipt_for_NB;  
        
        if nPT_from_SDi(NB)~=length(rdata1)
            error(' Received %1i points instead of the expected %1i pts from SD%1i',...
                  length(rdata1),nPT_from_SDi(NB),NB);
        else
            if verbose; fprintf(' Received the correct number of pts.\n'); end
        end
        if nPT_from_SDi(NB)
            j               = i + nPT_from_SDi(NB) - 1;
            els_D_NB(i:j)   = rdata1; % GLOBAL element number
            gX_PT_NB(:,i:j) = rdata2;
            ptr2nbPTs(i:j)  = NB;
            i               = j + 1;
        end
    else
        if verbose; fprintf(' Have to wait during this comm cycle...\n'); end
    end
end

end % END OF SUBFUNCTION communicate_pts

% #########################################################################

function els = guess_elems(EL2NOD)

[~,ind] = unique(EL2NOD);
els     = uint32(ceil(ind./size(EL2NOD,1)));

end % END OF SUBFUNCTION guess_elems