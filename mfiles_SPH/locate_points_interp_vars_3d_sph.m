function [vars_PT,els,locinfo] = locate_points_interp_vars_3d_sph...
    (MESH,COMM,gTH_PT,vars,OPTS_search,OPTS_interp,els)
% Usage: [vars_PT,els,locinfo] = locate_points_interp_vars_3d_sph...
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
% JH Feb 2017 : now works for any nmg (it only searches on the finest grid)

GCOORD     = MESH.GCOORD;
GCOORD_SPH = MESH.GCOORD_SPH;
EL2NOD     = MESH.EL2NOD;
nmg        = length(MESH.EL2NOD);
nnod       = size(GCOORD_SPH,2);
if size(vars,2)==nnod; vars=vars'; end
npt        = size(gTH_PT,2);
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
    check_lc = 1; % 1 --> calc local coords; find points to be located
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
    isrch           = find(isrch);
else
    isrch = 1:npt;
    lc    = zeros(3,npt);
    els   = ones(npt,1,'uint32');
end

tic
% Search points in coarse mesh
%     [els(isrch),gX_PT(:,isrch)] = tsearch3...
%         (GCOORD,EL2NOD{nmg}(1:4,:),gX_PT(:,isrch),OPTS_search,els(isrch));

%     [els(isrch),gTH_PT(:,isrch),OPTS_search] = tsearch2_sph_old_version...
%         (GCOORD_SPH,EL2NOD{1},gTH_PT(:,isrch),OPTS_search,MESH);
[els(isrch),gTH_PT(:,isrch),OPTS_search] = tsearch2_sph...
    (GCOORD_SPH,EL2NOD{1},gTH_PT(:,isrch),OPTS_search,MESH);
locinfo.t_srch       = toc;
locinfo.npt_srch     = length(isrch);
                                                                            tic
lc(:,isrch) = local_coords_3d_sph(GCOORD_SPH,EL2NOD{1},els(isrch),gTH_PT(:,isrch),MESH);
inot_loc    = any(lc(:,isrch)<0-tol_lc) | any(lc(:,isrch)>1+tol_lc) | ...
              sum(lc(:,isrch))>1+tol_lc;
inot_loc    = isrch(inot_loc);
if ~isempty(inot_loc)
    [els(inot_loc),lc(:,inot_loc)] = locate_points_in_neighboring_element...
        (els(inot_loc),gTH_PT(:,inot_loc),GCOORD_SPH,EL2NOD{1},MESH);
end

% Adjust local coordinates so that  0 <= [r,s,t,u=1-r-s-t] <= 1
lc(:,isrch)    = adjust_local_coords( lc(:,isrch) ); % *SUBFUNCTION*
locinfo.t_lc(1)=locinfo.t_lc(1)+toc;

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

% Interpolate variables at points
% ===============================
nnodel  = size(EL2NOD{1},1);
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
end

if strcmp(OPTS_interp.method_interp,'cubic')
    vars_PT = interp3d_cubic_p(GCOORD,EL2NOD{1},COMM,els,lc,vars,OPTS_interp); % interpolate
else
    vars_PT = interp3d_tet(EL2NOD{1},els,lc,vars,OPTS_interp); % interpolate
end

end % END OF FUNCTION locate_points_interp_vars_3d

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

function els = guess_elems(EL2NOD)

[~,ind] = unique(EL2NOD);
els     = uint32(ceil(ind./size(EL2NOD,1)));

end % END OF SUBFUNCTION guess_elems