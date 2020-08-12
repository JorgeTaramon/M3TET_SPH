function [MESH,COMM] = create_subdomain_multigrid_meshes_3d(MESH,COMM,nmg,fidl)
% Usage: [MESH,COMM] = create_subdomain_multigrid_meshes_3d(MESH,COMM,nmg,fidl)
%
% Purpose: Splits the domain into 'nsd' subdomains, creates all mutligrid
%          levels, calculates all data needed for communication on all
%          levels
%
% Input:
%   MESH : [structure] : all FE mesh data
%   COMM : [strcuture] : inter-subdomain communication data (not yet defined)
%   nmg  : [integer]   : number of multigrid levels
%   fidl : [integer]   : handle to log file
%
% Output:
%   MESH : [structure] : all FE mesh data
%   COMM : [strcuture] : inter-subdomain communication data (now complete)
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Oct 2012
% JH March 2013
% JH Feb 2016 : 
%

if numlabs==1
    FigNo_split   = 76;
    FigNo_balance = 67;
else
    FigNo_split   = 0;
    FigNo_balance = 0;
end

% =========================================================================
% DEFINE FUNCTIONS FOR ALL SUBDOMAIN COMMUNICATION
% =========================================================================
COMM = communication_function_handles(COMM); % *SUBFUNCTION*


% =========================================================================
% ASSIGN ELEMENTS TO EACH SUBDOMAIN
% =========================================================================
myid       = COMM.myid;
nsd        = COMM.nsd;
use_nesdis = 0; % switch for subdomain generation: 0 --> bi-section
                %                                  1 --> metis

% Assign elements to subdomains; For each element in MESH, el2sd 
% contains an integer value that defines the element's subdomain.
GCOORD_D  = MESH.GCOORD;
EL2NOD_D  = MESH.EL2NOD;
PointID_D = MESH.PointID;

% Calculate best split: el2sd will contains an integer value that defines
% each element's subdomain.
fprintf(fidl,'\n CREATING SUBDOMAINS...\n');
if use_nesdis
    % Use METIS to create subdomains
    fprintf(fidl,' Using METIS to calculate subdomains.\n');
    el2sd = el2sd_nesdis_3d(GCOORD_D,EL2NOD_D,nsd);
else
    % Use algorithm developed by Joerg Hasenclever
    fprintf(fidl,' Using bi-section algorithm to calculate subdomains.\n');

    % Define a percentage of shared nodes that you are willing to "pay"
    % for one less communication cycle. Use values between 5 and ~20.
    % INFO: A smaller value will try to reduce the number of shared 
    %       nodes on subdomain boundaries, but more connected 
    %       subdomains will be generated (more communication cycles).
    %       A higher value does the opposite: It tries to minimize the
    %       communication cycles but there qill be more shared nodes.
    fewer_NB_wght = 5;
    el2sd = el2sd_bisect_3d(GCOORD_D,EL2NOD_D,COMM,fewer_NB_wght,fidl);
end

% TESTING
% if isfield(MESH,'el2sd')
%     el2sd_mutils = MESH.el2sd;
%     nnod         = size(GCOORD_D,2);
%     allSDNB      = zeros(nsd-1,nsd);
%     for isd=1:nsd
%         NBs = el2sd_mutils( sum(ismember(EL2NOD_D(1:4,:)',EL2NOD_D(1:4,el2sd_mutils==isd)'),2)>0 );
%         NBs = unique(NBs(NBs~=isd));
%         nNB = length(NBs);
%         allSDNB(1:nNB,isd) = NBs(:);
%     end
% 
%     comm_scheme = pairwise_comm_scheme(allSDNB);
%     ncomm       = size(comm_scheme,1);
% 
%     nnodSD = zeros(nsd,1);
%     for isd=1:nsd
%         nnodSD(isd) = length( unique(EL2NOD_D(:,el2sd_mutils==isd) ) );
%     end
%     sumSDnods  = sum(nnodSD);
%     pct_shared = 100*(sumSDnods-nnod)/nnod;
%     fprintf(fidl,' DECOMP BY MUTILS     |  # shared nodes  | # comm cycles \n');
%     fprintf(fidl,'                      | %6i (%5.1f %%) |      %2i \n',...
%                      sumSDnods-nnod,pct_shared,ncomm);
% end


% =========================================================================
% SAVE DOMAIN MESH ON COAREST MULTIGRID LEVEL (needed for point search)
% =========================================================================
nnod_D         = max(EL2NOD_D(:));
MESH.nel_D     = size(EL2NOD_D,2);
MESH.nnod_D    = nnod_D;
MESH.GCOORD_D  = GCOORD_D(:,1:nnod_D);
MESH.EL2NOD_D  = EL2NOD_D;
MESH.PointID_D = MESH.PointID;
MESH.PhaseID_D = MESH.PhaseID;
MESH.el2sd     = el2sd;

if FigNo_split % PLOT SUBDOMAINS
    fprintf('\n Drawing subdomain configuration (might take a few sec)...');
    figure(FigNo_split);clf
%     if mod(nsd,2)==0
%         ind = [1:2:nsd;nsd-1:-2:1  ]; ind = ind(:)';
%     else
%         ind = [1:2:nsd;nsd-1:-2:1 0]; ind = ind(1:end-1);
%     end
    cmap     = jet(nsd);
%     cmap     = cmap(ind,:);
    TMP      = GCOORD_D;
    TMP(1,:) = TMP(1,:) - 5000;
    TMP(2,:) = TMP(2,:) - 5000;
    for isd=1:nsd
        simpplot(TMP',EL2NOD_D(:,el2sd==isd)','p(:,1)<-5000 & p(:,2)<-5000',cmap(isd,:),cmap(isd,:));
        hold on
    end
    TMP(2,:) = TMP(2,:) + 10000;
    for isd=1:nsd
        simpplot(TMP',EL2NOD_D(:,el2sd==isd)','p(:,1)<-5000 & p(:,2)>5000',cmap(isd,:),cmap(isd,:));
        hold on
    end
    TMP(1,:) = TMP(1,:) + 10000;
    for isd=1:nsd
        simpplot(TMP',EL2NOD_D(:,el2sd==isd)','p(:,1)>5000 & p(:,2)>5000',cmap(isd,:),cmap(isd,:));
        hold on
    end
    TMP(2,:) = TMP(2,:) - 10000;
    for isd=1:nsd
        simpplot(TMP',EL2NOD_D(:,el2sd==isd)','p(:,1)>5000 & p(:,2)<-5000',cmap(isd,:),cmap(isd,:));
        hold on
    end
    view(150,40); xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
    fprintf('done.\n');
    
%     % ALTERNATIVE PLOT
%     figure(FigNo_split+1);clf
%     cmap = lines(nsd);
%     subplot(1,4,1);
%     for isd=1:nsd
%         simpplot(GCOORD_D',EL2NOD_D(:,el2sd==isd)','p(:,1)<0 & p(:,2)<0',cmap(isd,:),cmap(isd,:));
%         hold on
%     end
%     view(130,24); xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
%     subplot(1,4,2);
%     for isd=1:nsd
%         simpplot(GCOORD_D',EL2NOD_D(:,el2sd==isd)','p(:,1)<0 & p(:,2)>0',cmap(isd,:),cmap(isd,:));
%         hold on
%     end
%     view(-60,24); xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
%     subplot(1,4,3);
%     for isd=1:nsd
%         simpplot(GCOORD_D',EL2NOD_D(:,el2sd==isd)','p(:,1)>0 & p(:,2)>0',cmap(isd,:),cmap(isd,:));
%         hold on
%     end
%     view(120,24); xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
%     subplot(1,4,4);
%     for isd=1:nsd
%         simpplot(GCOORD_D',EL2NOD_D(:,el2sd==isd)','p(:,1)>0 & p(:,2)<0',cmap(isd,:),cmap(isd,:));
%         hold on
%     end
%     view(-60,24); xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
    
%     % ALTERNATIVE PLOT
%     TMP.GCOORD     = GCOORD_D;
%     TMP.EL2NOD{1}  = EL2NOD_D;
%     TMP.nvertx     = 4;
%     TMP.PointID    = MESH.PointID;
%     TMP.DB_indices = MESH.DB_indices;
%     plot_domain_surfaces(FigNo_split,TMP,[1 2],el2sd,[]);
%     colormap(lines(nsd));
end


% =========================================================================
% EXTRACT THE SUBDOMAIN PART (el2sd==myid) FROM THE COARSEST MESH 
% =========================================================================
[GCOORD_SD,EL2NOD_SD,PointID_SD,nod_SD2D,nod_D2SD] = extract_SD_mesh...
   (GCOORD_D,EL2NOD_D,PointID_D,el2sd,myid,nsd); % *SUBFUNCTION*
MESH.EL2NOD                = cell(1,nmg);
MESH.EL2NOD{nmg}           = EL2NOD_SD;
MESH.nmg                   = nmg;
MESH.el_D2SDc              = zeros(1,MESH.nel_D,'uint32');
MESH.el_D2SDc(el2sd==myid) = 1:size(EL2NOD_SD,2);
if length(MESH.PhaseID)>1
    PhaseID_SD = MESH.PhaseID(el2sd==myid); % element PhaseID's
else
    PhaseID_SD = 1; %ones(1,MESH.nel,'int32') * MESH.PhaseID;
end

if FigNo_balance % PLOT LOAD BALANCE
    tmp = zeros(1,nsd);
    for isd=1:nsd
        tmp(isd) = length(nod_SD2D{isd});
    end
    figure(FigNo_balance);clf
    nnod_mean = mean(tmp);
    tmp       = 100 + (tmp-nnod_mean)./nnod_mean;
    subplot(2,1,1);bar(tmp);axis tight;xlabel('subdomain');ylabel('nnod');
    title(sprintf('average nnod per SD: %i',round(nnod_mean)));
    set(gca,'YLim',[min(tmp)-1 max(tmp)+1]);
    line(get(gca,'XLim'),[100 100],'Color','r','LineWidth',2);
    tmp      = accumarray(el2sd,ones(size(el2sd)));
    nel_mean = mean(tmp);
    tmp      = 100 + (tmp-nel_mean)./nel_mean;
    subplot(2,1,2);bar(tmp);axis tight;xlabel('subdomain');ylabel('nel');
    title(sprintf('average nel per SD: %i',round(nel_mean)));
    set(gca,'YLim',[min(tmp)-1 max(tmp)+1]);
    line(get(gca,'XLim'),[100 100],'Color','r','LineWidth',2);
end


% =========================================================================
% CALCULATE BOUNDARIES BETWEEN SUBDOMAINS
% =========================================================================
COMM = SD_boundaries_coarse_mesh(EL2NOD_SD,EL2NOD_D,el2sd,nod_SD2D,nod_D2SD,...
    COMM,nmg,fidl); % *SUBFUNCTION*


% =========================================================================
% RECURSIVELY REFINE SUBDOMAIN MESHES TO CREATE MULTIGRID LEVELS
% UPDATE INTER-SUBDOMAIN COMMUNICATION ARRAYS
% =========================================================================
ndof = 3; % number of degrees of freedom per node (velocity --> 3)
for img=nmg:-1:2
    % SAVE CURRENT MESH AS 'COARSE'
    EL2NOD_SDc  = EL2NOD_SD;
    GCOORD_SDc  = GCOORD_SD;
    PointID_SDc = PointID_SD;
    PhaseID_SDc = PhaseID_SD;
    
    % CALCULATE NEW REFINED MESH
    [GCOORD_SD,EL2NOD_SD,PointID_SD,PhaseID_SD,COMM] = generate_next_MG_level...
        (GCOORD_SDc,EL2NOD_SDc,PointID_SDc,PhaseID_SDc,MESH.DB_indices,COMM,img); % *SUBFUNCTION*
    
    % STORE CONNECTIVITY MATRIX
    MESH.EL2NOD{img-1} = EL2NOD_SD;
    
    % CALCULATE INTERMESH-TRANSFER OPERATOR (INTERPOLATION MATRIX)
    MESH.Ic2f{img-1}   = intermesh_transfer_matrix...
        (GCOORD_SD,EL2NOD_SD,GCOORD_SDc,EL2NOD_SDc,ndof); % *SUBFUNCTION*
end


% =========================================================================
% STORE ALL MESH DATA IN STRCUTURE 'MESH'
% CALCULATE TOTAL SIZE OF DOMAIN MESH (INFO ONLY) 
% =========================================================================
MESH.GCOORD  = GCOORD_SD;
MESH.nel     = size(MESH.EL2NOD{1},2);
MESH.nnod    = size(MESH.GCOORD,2);
MESH.PointID = PointID_SD;
MESH.PhaseID = PhaseID_SD;
for img=1:nmg
    COMM.nel_SD(img)  = size(MESH.EL2NOD{img},2);
    COMM.nnod_SD(img) = max(max(MESH.EL2NOD{img}));
    COMM.nel_D(img)   = MESH.nel_D * 8^(nmg-img);
    COMM.nnod_D(img)  = COMM.sum_all(length(COMM.unique_nodes{img}));
    COMM.shared_nodes{img} = ...
        setdiff(uint32(1:COMM.nnod_SD(img)),COMM.unique_nodes{img});
end
COMM.nod_SD2Dc = nod_SD2D;

end % END OF FUNCTION create_subdomain_multigrid_meshes_3d

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function [GCOORD_SD,EL2NOD_SD,PointID_SD,nod_SD2D,nod_D2SD] = ...
             extract_SD_mesh(GCOORD,EL2NOD,PointID,el2sd,myid,nsd)

nnodel   = size(EL2NOD,1);
nnod     = max(max(size(EL2NOD)));
nod_SD2D = cell(nsd,1);
nod_D2SD = cell(nsd,1);
for isd=1:nsd
    EL2NOD_SDi           = EL2NOD(:,el2sd==isd); % rows of "nodes" with elements of SD "isd"
    [nod_SDi2D,~,ind]    = unique(EL2NOD_SDi(:)); % unique list of the nodes in SD "isd"
    nel_SD               = size(EL2NOD_SDi,2);  % number of elements in SD "isd"
    nnod_SD              = length(nod_SDi2D); % number of nodes in SD "isd"
    nod_D2SDi            = zeros(nnod,1); % pointer global node numbers in domain
    nod_D2SDi(nod_SDi2D) = 1:nnod_SD;     % to global node numbers in subdomain "isd"
    % store domain to subdomain (global node numbers) pointer for SD "isd"
    nod_D2SD{isd}        = nod_D2SDi(:)';
    % store subdomain to domain (global node numbers) pointer for SD "isd"
    nod_SD2D{isd}        = nod_SDi2D(:)';

    if myid==isd
        % Store the subdomain mesh
        GCOORD_SD = GCOORD(:,nod_SDi2D);       % coordinates of all nodes in SD "isd
        EL2NOD_SD = reshape(uint32(ind),nnodel,nel_SD); % element connectivity in SD "isd"
        
        % Boundary nodes and domain boundary indices
        PointID_SD = PointID(nod_SDi2D);
    end
end

end % END OF SUBFUNCTION extract_SD_mesh

% #########################################################################

function COMM = SD_boundaries_coarse_mesh(EL2NOD_SD,EL2NOD_D,el2sd,nod_SD2D,nod_D2SD,COMM,nmg,fidl)

myid = COMM.myid;
nsd  = COMM.nsd;

% Calculate number of neighbors for all subdomains
allSDNB  = zeros(nsd-1,nsd);
for isd=1:nsd
    NBs = el2sd( sum(...
        ismember(EL2NOD_D(1:4,:),EL2NOD_D(1:4,el2sd==isd))...
                     ,1)>0 );
    NBs = unique(NBs(NBs~=isd));
    nNB = length(NBs);
    allSDNB(1:nNB,isd) = NBs(:);
end

% Calculate a pairwise communication scheme between all subdomains that
% share nodes. This pairwise scheme will avoids deadlocks during the
% subdomain-to-subdomain communication.
comm_scheme = pairwise_comm_scheme(allSDNB);
nNB         = size(comm_scheme,1);
NBlist      = comm_scheme(:,myid);

fprintf(fidl,'\n Communication scheme (%1i cycles): Subdomain...\n',nNB);
fprintf(fidl,repmat(' %2i',1,nsd),1:nsd);
fprintf(fidl,'\n talks to...\n');
for iNB=1:nNB
    frmt = [repmat(' %2i',1,nsd) '\n'];
    fprintf(fidl,frmt,comm_scheme(iNB,:));
end
fprintf(fidl,'\n');

% Allocate structure for communication data
COMM.nNB          = nNB;
COMM.NB           = NBlist';
COMM.comm_scheme  = comm_scheme;
COMM.mynod_SDB    = cell(1,nmg);
COMM.unique_nodes = cell(1,nmg);
for img=1:nmg
    COMM.mynod_SDB{img} = cell(1,nNB);
    COMM.nnod_SDB{img}  = zeros(1,nNB);
end

% Initialize a unique list of domain nodes (if all SDs merge there unique
% nodes, not a single domain node is missing and not a single one is
% counted twice; it's required, for instance, for parallel inner product)
nnod         = max(EL2NOD_SD(:));
unique_nodes = 1:nnod;
for iNB=1:nNB % Loop over subdomain neighbors
    NB = NBlist(iNB); % "name" or "index" of neighboring subdomain
    if NB==0
        % zero means that SD "myid" has to pause during this comm cycle
        continue
    end
    
    % Note: nod_SDB are the gobal node numbers in the domain (not the
    % subdomain node numbers)
    nod_SDB   = intersect(nod_SD2D{myid},nod_SD2D{NB});
    COMM.nnod_SDB{nmg}(iNB)  = length(nod_SDB);
    mynodSDB                 = nod_D2SD{myid}(nod_SDB);
    COMM.mynod_SDB{nmg}{iNB} = mynodSDB;
    
    % Update the unique node list (nodes on SD boundaries are assigned to
    % the subdomain with higher index)
    if myid<NB
        unique_nodes(mynodSDB) = 0;
    end
end
COMM.unique_nodes{nmg} = uint32(find(unique_nodes));

fprintf(fidl,'\n Communication scheme (%1i cycles)\n Subdomain...\n',nNB);
fprintf(fidl,repmat(' %2i',1,nsd),1:nsd);
fprintf(fidl,'\n talks to...\n');
for iNB=1:nNB
    frmt = [repmat(' %2i',1,nsd) '\n'];
    fprintf(fidl,frmt,comm_scheme(iNB,:));
end
fprintf(fidl,'\n');

end % END OF SUBFUNCTION SD_boundaries_coarse_mesh

% #########################################################################

function [GCOORD,EL2NOD,PointID,PhaseID,COMM] = generate_next_MG_level...
    (GCOORD_c,EL2NOD_c,PointID_c,PhaseID_c,DB_indices,COMM,img)

FigNo  = 0;  % set to zero to NOT show mesh figures

if FigNo
    figure(FigNo);clf;
    tetramesh(EL2NOD_c(1:4,:)',GCOORD_c',...
              'EdgeColor','b',...
              'FaceAlpha',0.1,'FaceColor','k',...
              'LineStyle','-','LineWidth',1);
end

% =========================================================================
% Check whether the mesh has linear (4-node) or quadratic (10-node) 
% tetrahedral elements. If quadratic, split each element into 8 linear
% sub-elements.
% =========================================================================
nnodel = size(EL2NOD_c,1);
if nnodel==10
    [EL2NOD_c,PhaseID_c] = tetmesh_p2_to_p1(EL2NOD_c,PhaseID_c);
end

% =========================================================================
% Convert 4-node tetrahedra mesh into 10-node mesh by creating nodes on 
% each tetrahedron's edges.
% =========================================================================
[GCOORD,EL2NOD,PointID] = tetmesh_p1_to_p2(GCOORD_c,EL2NOD_c,PointID_c,DB_indices);

% =========================================================================
% Update the communication arrays for the newly created nodes
% =========================================================================
if numlabs>1
    COMM = calc_COMM_p2(GCOORD,EL2NOD,COMM,img);
end

% =========================================================================
% If the coarse mesh had linear 4-node tetrahedra, re-connect nodes the
% current 10-node tetra mesh to obtain a 4-node mesh again.
% =========================================================================
if nnodel==4
    [EL2NOD,PhaseID] = tetmesh_p2_to_p1(EL2NOD,PhaseID_c);
else
    PhaseID          = PhaseID_c;
end

end % END OF SUBFUNCTION generate_next_MG_level

% #########################################################################

function COMM = calc_COMM_p2(GCOORD,EL2NOD,COMM,img)

save_mesh_to_tecfile = 0;

% =========================================================================
% Get boundary index information for all new nodes
% =========================================================================
nnod  = max(EL2NOD(:));

% Create a pointer that for each edge node in the mesh returns the two
% bounding vertex nodes
pointer_Enod2Vnod                 = zeros(nnod,2,'uint32');
pointer_Enod2Vnod(EL2NOD( 5,:),:) = EL2NOD([1 2],:)';
pointer_Enod2Vnod(EL2NOD( 6,:),:) = EL2NOD([2 3],:)';
pointer_Enod2Vnod(EL2NOD( 7,:),:) = EL2NOD([3 4],:)';
pointer_Enod2Vnod(EL2NOD( 8,:),:) = EL2NOD([1 4],:)';
pointer_Enod2Vnod(EL2NOD( 9,:),:) = EL2NOD([1 3],:)';
pointer_Enod2Vnod(EL2NOD(10,:),:) = EL2NOD([2 4],:)';
pointer_Enod2Vnod                 = sort(pointer_Enod2Vnod,2);

% Initialize list of unique nodes
unique_nodes = uint32(1:nnod);
for iNB=1:COMM.nNB
    NB        = COMM.NB(iNB);
    if NB==0
        continue
    end
    
    mynod_SDB_p1 = COMM.mynod_SDB{img}{iNB};
    
    % CHECK SUBDOMAIN BOUNDARY
    data_to_NB   = single(GCOORD(:,mynod_SDB_p1)');
    data_from_NB = COMM.sendnrecv_1NB(NB,data_to_NB);
    if max(max(abs(data_to_NB-data_from_NB)))>1e-8
        error(' This must not happen');
%     else
%         disp(' SUBDOMAIN BOUNDARY NODES ARE IN CORRRECT ORDER.');
    end
    % CHECK SUBDOMAIN BOUNDARY

    % =====================================================================
    % Assume first that each new edge node that is bounded by two vertex 
    % nodes that are already shared with the neighboring subdomain is also 
    % shared with this neighbor. Note that this is not always true (see 
    % the example below) and a second check is required afterwards!
    % =====================================================================
    %
    %           c
    %          / \
    %         /   \        subdomain 1 (SD1), FACE of tetrahedra 1
    %        /  1  \
    %       /       \
    %      a----x----b
    %      A         B     subdomain 3 (SD3), FACE of tetrahedra 2
    %     / \       / \
    %    /   \  2  /   \   subdomain 2 (SD2), FACES of tetrahedra 3 and 4
    %   /     \   /     \
    %  /   3   \d/   4   \
    % E---------D---------F
    % 
    % For instance, if nodes a,b,d of SD1 (top, lower case letters) are
    % shared with nodes A,B,D of SD2 (bottom, capital letters). The new
    % edge node "x" in SD1 won't necessarily be shared with SD2, even 
    % though it is located between two shared nodes (a, b). The element #2
    % could belong to a 3rd subdomain... Only SD2 knows whether or not it 
    % also contains a node in thelocation of node "x" (SD2 knows if it
    % created a new edge node between A and B). This will be the second
    % check now.
    
    % Find edge nodes that are bounded by shared vertex nodes:
    [tf,~]       = ismember(pointer_Enod2Vnod,mynod_SDB_p1);
    % The edge nodes possibly shared with the neighbor are:
    mynod_SDB_p2 = find(sum(tf,2)==2); % result will always be sorted!!!
    
    % Check 2: Compare the coordinates of the edge nodes of which NB thinks
    % they are shared
    if numlabs>1
        GCOORD_SD = single(GCOORD(:,mynod_SDB_p2)');
        if COMM.myid<NB
            COMM.send_1NB(NB,GCOORD_SD,333);
            loc_SD            = COMM.recv_1NB(NB,334);
        else
            GCOORD_NB         = COMM.recv_1NB(NB,333);
            [~,loc_SD,loc_NB] = intersect(GCOORD_SD,GCOORD_NB,'rows');
            COMM.send_1NB(NB,loc_NB,334);
        end
        mynod_SDB_p2 = mynod_SDB_p2(loc_SD);
    end
% % % %     % Now send the neighbor the vertex nodes that bound all edge nodes that
% % % %     % are possible shared. The location of these vertex nodes in list
% % % %     % "mynod_SDB_p1" will be communicated.
% % % %     loc_Vnods_bounding_shared_Enods = loc(mynod_SDB_p2,:);
% % % %     
% % % %     if numlabs>1
% % % %         data_from_NB = COMM.sendnrecv_1NB(NB,loc_Vnods_bounding_shared_Enods);
% % % %     else
% % % %         data_from_NB = loc_Vnods_bounding_shared_Enods;
% % % %     end
    
% %     if numlabs>1
% %         filename = ['Workspace1_Lab' num2str(labindex)];
% %         save(filename);
% %     else
% %         myid     = 2;
% %         filename = ['Workspace1_Lab' num2str(myid)];
% %         load(filename);
% %     end

% % % %     % Form a list of pairs of vertex nodes between the NB thinks it shares
% % % %     % an edge node. 
% % % %     Vnods_bounding_shared_Enods = mynod_SDB_p1( data_from_NB );
% % % %     Vnods_bounding_shared_Enods = sort(Vnods_bounding_shared_Enods,2);
% % % %     
% % % %     % CHECK 2: Was an edge node created in the location where the neighbor
% % % %     % thinks it shares an edge node? I.e. (see above) did SD 2 create a
% % % %     % node between A,B where SD 1 thinks it shares edge node x because it
% % % %     % is located between already share vertex nodes a,b?
% % % %     [tf,mynod_SDB_p2] = ismember(Vnods_bounding_shared_Enods,pointer_Enod2Vnod,'rows');
% % % %     
% % % %     % Nodes that are not shared have to be communicated. They will be
% % % %     % deleted from the list of subdomain boundary nodes. Both subdomains
% % % %     % have to do this.
% % % %     % Index of nodes where SD knows that they are not shared:
% % % %     loc_not_shared1 = find(~tf)
% % % %     % Send "loc_not_shared1" and receive index of nodes where NB knows that
% % % %     % they are not shared:
% % % %     loc_not_shared2 = COMM.sendnrecv_1NB(NB,loc_not_shared1) 
% % % %     
% % % %     if numlabs>1
% % % %         filename = ['Workspace2_Lab' num2str(labindex)];
% % % %         save(filename);
% % % %     else
% % % %         myid     = COMM.myid;
% % % %         filename = ['Workspace2_Lab' num2str(myid)];
% % % %         load(filename);
% % % %     end
% % % %     
% % % %     loc_not_shared = union(loc_not_shared1,loc_not_shared2)
% % % %     mynod_SDB_p2(loc_not_shared) = [];
% % % %     length(mynod_SDB_p2)
% % % %     
% % % %     stop
% % % %     
% % % %     % Now the new edge nodes shared with neighbor NB are known but they
% % % %     % still need to be sorted.
% % % %     if numlabs>1
% % % %         if COMM.myid<NB
% % % %             mynod_SDB_p2 = sort(mynod_SDB_p2);
% % % %             data_to_NB   = single(GCOORD(:,mynod_SDB_p2)');
% % % %             COMM.send_1NB(NB,data_to_NB,333);
% % % %         
% % % %             filename = ['Workspace3_Lab' num2str(labindex)];
% % % %             save(filename);
% % % %             stop
% % % %         else
% % % %             data_from_NB = COMM.recv_1NB(NB,333);
% % % %             filename = ['Workspace3_Lab' num2str(labindex)];
% % % %             save(filename);
% % % %             stop
% % % %             
% % % %             [~,perm]     = ismember(single(GCOORD(:,mynod_SDB_p2)'),data_from_NB,'rows');
% % % %             mynod_SDB_p2 = mynod_SDB_p2(perm);
% % % %         end
% % % %     else
% % % %         myid     = 2;
% % % %         if COMM.myid<NB
% % % %             filename = ['Workspace3_Lab' num2str(myid)];
% % % %             load(filename);
% % % %         else
% % % %             filename = ['Workspace3_Lab' num2str(myid)];
% % % %             load(filename);
% % % %             
% % % %             [~,perm]     = ismember(single(GCOORD(:,mynod_SDB_p2)'),data_from_NB,'rows');
% % % %             mynod_SDB_p2 = mynod_SDB_p2(perm);
% % % %         end
% % % %     end
% % % %     
% % % %     stop
    
    % CHECK SUBDOMAIN BOUNDARY
    data_to_NB   = single(GCOORD(:,mynod_SDB_p2)');
    data_from_NB = COMM.sendnrecv_1NB(NB,data_to_NB);
    if max(max(abs(data_to_NB-data_from_NB)))>1e-8
        error(' This must not happen');
%     else
%         disp(' SUBDOMAIN BOUNDARY NODES ARE IN CORRRECT ORDER.');
    end
    % CHECK SUBDOMAIN BOUNDARY
    
    mynod_SDB_p2               = [mynod_SDB_p1(:)' mynod_SDB_p2(:)'];
    COMM.mynod_SDB{img-1}{iNB} = mynod_SDB_p2;
    COMM.nnod_SDB{img-1}(iNB)  = length(mynod_SDB_p2);
    
    if COMM.myid<NB
        unique_nodes( mynod_SDB_p2 ) = 0;
    end
end

COMM.unique_nodes{img-1} = uint32(find(unique_nodes));

% if numlabs>1
%     filename = ['Workspace5_Lab' num2str(labindex)];
%     save(filename);
% else
%     myid     = 1;
%     filename = ['Workspace5_Lab' num2str(myid)];
%     load(filename);
% end


% For debugging: write subdomain boundaruy data to Tecplot file 
% =============================================================
if save_mesh_to_tecfile
    tecfile = ['MESH_SD' num2str_d(COMM.myid,2) '.dat'];
    nnod = max(EL2NOD(:));
    s    = cell(1,4);
    s{1} = ['VARIABLES = "X", "Y", "Z"' sprintf(', "NB%1i"',1:COMM.nNB) '\n'];
    s{2} = sprintf('ZONE T="%s",N=%1i, E=%1i, DataPacking=POINT, ZoneType=FETETRAHEDRON \n',...
              ['SD' num2str_d(COMM.myid,2)],nnod,8*size(EL2NOD,2));
    data = zeros(nnod,COMM.nNB);
    for iNB=1:COMM.nNB
        NB = COMM.NB(iNB);
        if NB==0
            continue
        end
        mynod_SDB           = COMM.mynod_SDB{img-1}{iNB};
        data(mynod_SDB,iNB) = NB;
    end
    data = [GCOORD' data]';
    fmt  = repmat('%6E ',1,3+COMM.nNB); fmt(end:end+1) = '\n';
    s{3} = sprintf(fmt,data);

    s{4} = sprintf('%1i %1i %1i %1i\n',double(tetmesh_p2_to_p1(EL2NOD)));
    output_unit = fopen(tecfile,'w');
    fprintf(output_unit, [s{1:4}]);
    fclose(output_unit);
end

end % END OF SUBFUNCTION calc_COMM_p2

% #########################################################################

function [GCOORD_SD,EL2NOD_SD,DBnods_SD,COMM,el2face,edge2sd,face2sd] = ...
    generate_next_MG_level_old(GCOORD_SDc,EL2NOD_SDc,DBnods_SDc,DB_indices,...
                               COMM,img,el2face_c,edge2sd_c,face2sd_c)

myid = COMM.myid;
nsd  = COMM.nsd;

% size of coarse mesh
nel_c      = size(EL2NOD_SDc,2); % number of elements
nnod_c     = size(GCOORD_SDc,2); % number of nodes
EL2NOD_SDc = EL2NOD_SDc'; % connectivity matrix

FigNo = 0; % img*1000 + myid*100;  % set to zero to NOT show mesh figures
if FigNo>0
    view3D = [20 20];
    switch myid
        case 1
            view3D = [-12,12];
            view3D = [-146 26];
            view3D = [20 20];
        case 2
            view3D = [-76,12];
            view3D = [-12,12];
            view3D = [-73 -24];
        case 3
            view3D = [70 10];
        case 4
            view3D = [-145 22];
    end
    fs    = 0;
    flag  = 4;
    fc    = 'k';
    fa    = 0.8;
    lc    = 'b';
    lw    = 1;
%     sfigure(FigNo);clf
%     showmesh_3D(FigNo,sprintf('MG%1i SD%1i of %1i ',img,myid,nsd),...
%                 GCOORD_SDc,EL2NOD_SDc,fs,flag,fc,fa,lc,lw);
%     view(view3D);
%     saveas(gcf,sprintf('%s/Mesh_MG%1i_SD%1ix%1i',outdir,img,myid,nsd),'fig');
end
% open ../Output/Test_p/Mesh_MG2_SD1x2.fig
% open ../Output/Test_p/Mesh_MG2_SD2x2.fig

% (1) Create bars that connect nodes between which a new node needs to be 
%     generated
% =======================================================================

% Create a pointer bars defined by their end-nodes
% (e.g. bar 1 has end-nodes [1 5], bar 2 has [5 2],...)
bar2node = reshape(EL2NOD_SDc(:,[ 1  5  5  2 ... % 1st edge of parent
                                2  6  6  3 ... % 2nd edge of parent
                                3  7  7  4 ... % 3rd edge of parent
                                4  8  8  1 ... % 4th edge of parent
                                1  9  9  3 ... % 5th edge of parent
                                4 10 10  2 ... % 6th edge of parent
                               10  6  6  7  7 10 ... % 3 edges on face #1 of parent (face opposing node 1)
                                9  8  8  7  7  9 ... % 3 edges on face #2 of parent (face opposing node 2)
                                8  5  5 10 10  8 ... % 3 edges on face #3 of parent (face opposing node 3)
                                5  9  9  6  6  5 ... % 3 edges on face #4 of parent (face opposing node 4)
                                6  8 ...       % node in center of parent
                                ])',2,[] )';

% Find the bars that are shared by neighboring elements and return a unique
% list of bars. The bars are the new edges in the refined mesh!!!!
[bar2node,~,ib] = unique_keep_order(bar2node);
el2bar = reshape(ib,25,nel_c)'; % element to bar connectivity after doubles 
                                % have been merged (removed)
nbars  = size(bar2node,1); % number of unique bars (i.e. number of new edge 
                           % nodes that have to be generated)

% Coordinates of the nodes defining each bar's end points
xBarEnds = reshape(GCOORD_SDc(1,bar2node'),2,[]);
yBarEnds = reshape(GCOORD_SDc(2,bar2node'),2,[]);
zBarEnds = reshape(GCOORD_SDc(3,bar2node'),2,[]);

% Create new node at each bar mid-point
xBarMids = 0.5*sum(xBarEnds,1);
yBarMids = 0.5*sum(yBarEnds,1);
zBarMids = 0.5*sum(zBarEnds,1);
iBarMids = nnod_c + (1:nbars); % global node number of new nodes

% (2) Allocate storage for mesh
% =============================
nel    = 8*nel_c;
nnod   = nnod_c + nbars;
nVnod  = nnod_c;
nEnod  = nbars;
GCOORD_SD                  = zeros(3,nnod); % storage for fine mesh node coordinates
GCOORD_SD(:,1:nnod_c)      = GCOORD_SDc;
GCOORD_SD(1,nnod_c+1:nnod) = xBarMids;
GCOORD_SD(2,nnod_c+1:nnod) = yBarMids;
GCOORD_SD(3,nnod_c+1:nnod) = zBarMids;

% if FigNo
%     for iel_c=1:nel_c
%         sfigure(FigNo);clf;
%         patch_tetrahedron(FigNo,GCOORD_SDc,EL2NOD_SDc,iel_c,3,'r','-',1,fontsize);
%     end
% end

% (3) New connectivity matrix
% (3a) Split each coarse mesh 10-node element into 8 linear sub-elements
%      (all nodes on the coarse mesh become vertex nodes of the fine mesh)
% I.e.: element #1 has nodes 1 5 9  8
%       element #2 has nodes 5 2 6 10
%       etc
EL2NOD_SD        = zeros(nel,10); % storage for fine mesh connectivity matrix
EL2NOD_SD(:,1:4) = reshape(EL2NOD_SDc(:,[ 1  5  9  8 ...
                                     5  2  6 10 ...
                                     9  6  3  7 ...
                                     8 10  7  4 ...
                                     8  5  9 10 ...
                                     7  6 10  9 ...
                                     5  9 10  6 ...
                                     8 10  9  7 ])',4,[])';

% (3b) Write the edge nodes into the connectivity matrix
% Part 1: all nodes on edges on parent element
ibar  = [1 2 3 4 5 6 7 8 9 10 11 12]; % local number of bar (1:24)
subel = [1 2 2 3 3 4 4 1 1  3  4  2]; % sub-element in each parent element (1:8)
nod_f = [5 5 6 6 7 7 8 8 9  9 10 10]; % local node in sub-element (1:10)

% Part 2: all nodes on face #1 of parent element
ibar  = [ibar  13 13 13 14 14 15 15 15]; % local number of bar (1:24)
subel = [subel  2  6  7  3  6  4  6  8]; % sub-element in each parent element (1:8)
nod_f = [nod_f  7  6  7 10  5  6  9 10]; % local node in sub-element (1:10)

% Part 3: all nodes on face #2 of parent element
ibar  = [ibar  16 16 16 17 17 18 18 18]; % local number of bar (1:24)
subel = [subel  1  5  8  4  8  3  6  8]; % sub-element in each parent element (1:8)
nod_f = [nod_f  7  9  9  9  8  8  8  7]; % local node in sub-element (1:10)
% ibar  = [ibar  16 16 16 17 17 17 18 18]; % local number of bar (1:24)
% subel = [subel  1  5  8  8  3  6  4  8]; % sub-element in each parent element (1:8)
% nod_f = [nod_f  7  9  9  8  8  8  9  7]; % local node in sub-element (1:10)

% Part 4: all nodes on face #3 of parent element
ibar  = [ibar  19 19 20 20 20 21 21 21]; % local number of bar (1:24)
subel = [subel  1  5  2  5  7  4  5  8]; % sub-element in each parent element (1:8)
nod_f = [nod_f 10  5  8 10  9  5  8  5]; % local node in sub-element (1:10)

% Part 5: all nodes on face #4 of parent element
ibar  = [ibar  22 22 22 23 23 23 24 24]; % local number of bar (1:24)
subel = [subel  1  5  7  3  6  7  2  7]; % sub-element in each parent element (1:8)
nod_f = [nod_f  6  6  5  5 10 10  9  8]; % local node in sub-element (1:10)

% Part 6: node at center of parent element
ibar  = [ibar  25 25 25 25]; % local number of bar (1:24)
subel = [subel  5  6  7  8]; % sub-element in each parent element (1:8)
nod_f = [nod_f  7  7  6  6]; % local node in sub-element (1:10)

% Use the above triplet to write the connectivity for the edge nodes of all
% elements in the fine mesh
EL2NOD_SD = write_connect(EL2NOD_SD,ibar,subel,nod_f);
function EL2NOD_SD = write_connect(EL2NOD_SD,ibar,subel,nod_f)
    for jb=1:length(ibar)
        % g_bar is global bar number
        g_bar                   = el2bar(:,ibar(jb));
        % use subelement in each coarse mesh element (i.e. 1,2,..., or 8)
        % to calculate global number of element in fine mesh
        els_f                   = 8*(1:nel_c)'-8 + subel(jb)*ones(nel_c,1);
        % New global mid-side nodes in the fine mesh are given by iBarMids
        % (see above: iBarMids = nnod_c + (1:nbars);)
        EL2NOD_SD(els_f,nod_f(jb)) = iBarMids(g_bar);
    end
end

% if FigNo
%     for iel_c=1
%         sfigure(FigNo-2);clf;
%         patch_tetrahedron(FigNo-2,GCOORD_SDc,EL2NOD_SDc,iel_c,0,'r','-',1,fontsize);
%         for ii=1:25
%             text(xBarMids(el2bar(iel_c),ii),...
%                  yBarMids(el2bar(iel_c),ii),...
%                  zBarMids(el2bar(iel_c),ii),...
%                  sprintf('%1i',ii),'Fontsize',8);
%         end
%         sfigure(FigNo-1);clf;
%         patch_tetrahedron(FigNo-1,GCOORD_SDc,EL2NOD_SDc,iel_c,4,'r','-',1,fontsize);
%         sfigure(FigNo);clf;
%         patch_tetrahedron(FigNo,GCOORD_SDc,EL2NOD_SDc,iel_c,3,'r','-',1,fontsize);
%         iel_f = 8*(iel_c-1);
%         for ii=1:8
%             iel_f = iel_f+1;
%             sfigure(FigNo+ii);clf;
%             patch_tetrahedron(FigNo+ii,GCOORD_SDc,EL2NOD_SDc,iel_c,0,'k','--',0.5,fontsize)
%             patch_tetrahedron(FigNo+ii,GCOORD_SD,EL2NOD_SD,iel_f,3,'r','-',1,fontsize)
%             %title(sprintf('Element %1i (sub-element %1i)',iel_f,ii));
%             sfigure(2*FigNo+ii);clf;
%             patch_tetrahedron(2*FigNo+ii,GCOORD_SDc,EL2NOD_SDc,iel_c,0,'k','--',0.5,fontsize)
%             patch_tetrahedron(2*FigNo+ii,GCOORD_SD,EL2NOD_SD,iel_f,4,'r','-',1,fontsize)
%             %title(sprintf('Element %1i (sub-element %1i)',iel_f,ii));
%         end
%         1;
%     end
% end

% (5) Get boundary index information for new nodes
% ================================================
% (5a) Create a pointer that for each edge node in the mesh returns the 2
%      bounding vertex nodes
pointer_Edge2Vnod = zeros(nEnod,2);
pointer_Edge2Vnod(EL2NOD_SD(:, 5)-nVnod,:) = EL2NOD_SD(:,[1 2]);
pointer_Edge2Vnod(EL2NOD_SD(:, 6)-nVnod,:) = EL2NOD_SD(:,[2 3]);
pointer_Edge2Vnod(EL2NOD_SD(:, 7)-nVnod,:) = EL2NOD_SD(:,[3 4]);
pointer_Edge2Vnod(EL2NOD_SD(:, 8)-nVnod,:) = EL2NOD_SD(:,[4 1]);
pointer_Edge2Vnod(EL2NOD_SD(:, 9)-nVnod,:) = EL2NOD_SD(:,[1 3]);
pointer_Edge2Vnod(EL2NOD_SD(:,10)-nVnod,:) = EL2NOD_SD(:,[2 4]);


% (5b) Find fine mesh vertex nodes (i.e. coarse mesh nodes) that are in the
%      list of coarse domain boundary nodes (1st column in DBnods_SDc)
[tf,loc] = ismember(pointer_Edge2Vnod,DBnods_SDc(:,1)); % find DBnod_c in pointer_Edge2Vnod
tf       = find(tf);
% Create new pointer to the domain boundary index (2nd column in DBnods_SDc)
pointer_Edge2DBnod     = zeros(nEnod,2); % create new pointer
pointer_Edge2DBnod(tf) = DBnods_SDc(loc(tf),2);

% (5c) Use the pointer pointer_Edge2DBnod to check (for each NEW node) 
%      what the domain boundary index of its two neighborung vertex nodes
%      are. The domain boundary index for the new nodes can be found using
%      the following logic:
%      Edge nodes cannot be generated in domain corners but only on
%      domain edges or faces.
%      Face nodes cannot be generated on domain edges (and of course not in
%      domain corners).
%
%    108---211---107    % Coordinate system to define xmin, xmax, ymin, ymax, zmin, zmax
%     /:         /|     %     y
%  212 :  306  210|     %   /
%   /  :       /  |     %  /
%105---:209---106 |     % 0------- x
%  |  208  304|  207    % |
%  |   :      |   |     % |
%  |305:      |303|     % |
% 205  : 302 206  |     %-z
%  | 104---203|---103
%  |  /       |  /
%  |204  301  |202
%  |/         |/
%101--- 201---102

% Faces: 301 bottom
%        302 front
%        303 right
%        304 back
%        305 left
%        306 top
% Edges: 201 lower front
%        202 lower right
%        203 lower back
%        204 lower left
%        205 front left
%        206 front right
%        207 back right
%        208 back left
%        209 top front
%        210 top right
%        211 top back
%        212 top left
% Corners: 101 lower front left
%          102 lower front right
%          103 lower back right
%          104 lower back left
%          105 top front left
%          106 top front right
%          107 top back right
%          108 top back left
%
% special: 310 face nodes west of MAR
%          220 edge nodes west of MAR
%          311 face nodes east of MAR
%          221 edge nodes east of MAR
%          222 nodes on MAR
%          223 nodes in TF
%          110 nodes at JCT (intersection MAR-TF)

% There are 9 standard cases and 6 special cases in 3D:
%       1) both on same F ==> F
%       2) both on same E ==> E
%       3) C and E ==> E (C+E must be on same domain face)
%       4) C and F ==> F (C+F must be on same domain face)
%       5) E and F ==> F (E+F must be on same domain face)
%       6) C and 0 ==> 0
%       7) E and 0 ==> 0
%       8) F and 0 ==> 0
%       9) 0 and 0 ==> 0 
%      S1) C1 and C2 ==> E between C1 and C2
%      S2) E1 and E2 ==> F between E1 and E2
%      S3) F1 and F2 ==> 0 (inside domain)
%      S4) C1 and E2 ==> 0 (inside domain)
%      S5) C1 and F2 ==> 0 (inside domain)
%      S6) E1 and F2 ==> 0 (inside domain)
%      The numbers in the special cases indicate DIFFERENT domain faces !!!

% % % For testing and understanding the logic:
% % pointer_Edge2DBnod = [101 101; % case 1
% %                       201 201; % case 2
% %                       301 201; % case 3
% %                       301 101; % case 4
% %                       201 101; % case 5
% %                       301  0 ; % case 6
% %                       201  0 ; % case 7
% %                       101  0 ; % case 8
% %                        0   0 ; % case 9
% %                       301 302; % special case 1
% %                       201 202; % special case 2
% %                       101 102];% special case 3
% %                       301 202; % special case 4
% %                       301 103; % special case 5
% %                       201 103];% special case 6

% Cases 1:9 are taken care of by chosing the LARGER index of the 2
% bounding vertex nodes. Treat all nodes as standard cases first:
DBindx  = max( pointer_Edge2DBnod,[],2 );
% HOWEVER, if any index is zero, the node will be inside the domain. Set
% these indices euqla to zero now.
ind_0   = any(pointer_Edge2DBnod==0,2);
DBindx(ind_0) = 0;

% Need special treatment for the special cases...
% (a) all potential special cases
indS  = find( pointer_Edge2DBnod(:,1)~=pointer_Edge2DBnod(:,2) & ...
              pointer_Edge2DBnod(:,1)>0 & pointer_Edge2DBnod(:,2)>0 );
% (b) special case 1 (between 2 corners); all possible combinations
indS1 = indS(all(pointer_Edge2DBnod(indS,:)>100,2) & all(pointer_Edge2DBnod(indS,:)<200,2));
if ~isempty(indS1)
    ind_0          = indS1;
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[101 102]) ,2) );
    DBindx(indS1a) = 201; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[102 103]) ,2) );
    DBindx(indS1a) = 202; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[103 104]) ,2) );
    DBindx(indS1a) = 203; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[101 104]) ,2) );
    DBindx(indS1a) = 204; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[101 105]) ,2) );
    DBindx(indS1a) = 205; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[102 106]) ,2) );
    DBindx(indS1a) = 206; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[103 107]) ,2) );
    DBindx(indS1a) = 207; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[104 108]) ,2) );
    DBindx(indS1a) = 208; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[105 106]) ,2) );
    DBindx(indS1a) = 209; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[106 107]) ,2) );
    DBindx(indS1a) = 210; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[107 108]) ,2) );
    DBindx(indS1a) = 211; ind_0 = setdiff(ind_0,indS1a);
    indS1a         = indS1( all( ismember(pointer_Edge2DBnod(indS1,:),[105 108]) ,2) );
    DBindx(indS1a) = 212; ind_0 = setdiff(ind_0,indS1a);
    DBindx(ind_0)  = 0;
end

% (c) special case 2 (between 2 different edges); all possible combinations
indS2 = indS(all(pointer_Edge2DBnod(indS,:)>200,2) & all(pointer_Edge2DBnod(indS,:)<300,2));
if ~isempty(indS2)
    ind_0          = indS2;
    indS2a         = indS2( all( ismember(pointer_Edge2DBnod(indS2,:),[201 202 203 204]) ,2) );
    DBindx(indS2a) = 301; ind_0 = setdiff(ind_0,indS2a);
    indS2a         = indS2( all( ismember(pointer_Edge2DBnod(indS2,:),[201 205 206 209]) ,2) );
    DBindx(indS2a) = 302; ind_0 = setdiff(ind_0,indS2a);
    indS2a         = indS2( all( ismember(pointer_Edge2DBnod(indS2,:),[202 206 207 210]) ,2) );
    DBindx(indS2a) = 303; ind_0 = setdiff(ind_0,indS2a);
    indS2a         = indS2( all( ismember(pointer_Edge2DBnod(indS2,:),[203 207 208 211]) ,2) );
    DBindx(indS2a) = 304; ind_0 = setdiff(ind_0,indS2a);
    indS2a         = indS2( all( ismember(pointer_Edge2DBnod(indS2,:),[204 205 208 212]) ,2) );
    DBindx(indS2a) = 305; ind_0 = setdiff(ind_0,indS2a);
    indS2a         = indS2( all( ismember(pointer_Edge2DBnod(indS2,:),[209 210 211 212]) ,2) );
    DBindx(indS2a) = 306; ind_0 = setdiff(ind_0,indS2a);
    DBindx(ind_0)  = 0;
end

% (d) special case 3
indS3 = indS(all(pointer_Edge2DBnod(indS,:)>300,2) & all(pointer_Edge2DBnod(indS,:)<400,2));

% (e) special cases 4, 5 and 6
indS456 = setdiff(indS,[indS1(:);indS2(:);indS3(:)]);
% Now loop over the domain faces and check if the different DB indices
% (e.g. C1 and E2) are on the same domain face as defined by the cell-matrix 
% "DB_indices". If they are on teh same face, they are NOT INSIDE the
% domain!
for iDface=1:length(DB_indices)
    iSameFace = sum(ismember( pointer_Edge2DBnod(indS456,:),DB_indices{iDface} ),2)==2;
    indS456(iSameFace) = [];
end
% old version:
% Dfaces  = [301 201 202 203 204 101 102 103 104; % domain bottom     (zmin)
%            302 201 205 206 209 101 102 105 106; % domain front      (ymin)
%            303 202 206 207 210 102 103 106 107; % domain right side (xmax)
%            304 203 207 208 211 103 104 107 108; % domain back       (ymax)
%            305 204 205 208 212 101 104 105 108; % domain left side  (xmin)
%            306 209 210 211 212 105 106 107 108];% domain top        (zmax)
% % Now loop over the domain faces and check if the different DB indices
% % (e.g. C1 and E2) are on the same domain face as defined by matrix Dfaces.
% for iDface=1:size(Dfaces,1)
%     iSameFace = sum(ismember( pointer_Edge2DBnod(indS456,:),Dfaces(iDface,:) ),2)==2;
%     indS456(iSameFace) = [];
% end

% Special case 3,4,5 and 6 (between 2 different domain faces thus INSIDE
% the domain) ===> simply remove these nodes from the list of boundary nodes
inod_DB = (nVnod+1:nnod)';
DBindx ([indS3(:);indS456(:)]) = [];
inod_DB([indS3(:);indS456(:)]) = [];

% Remove all nodes that have a zero DB index
ind_0          = ~DBindx;
DBindx(ind_0)  = [];
inod_DB(ind_0) = [];

% merge coarse-mesh and new domain boundary nodes
DBnods_SD      = [DBnods_SDc; [inod_DB DBindx]];


if FigNo>0
    FigNo = FigNo + 1;
    sfigure(FigNo);clf;if numlabs>1;set(gcf,'visible','off'); end
    showmesh_3D(FigNo,sprintf('MG%1i SD%1i of %1i ',img+1,myid,nsd),...
                GCOORD_SD,EL2NOD_SD,fs,flag,fc,fa,lc,lw);
    view(view3D);
end

% if img==2
%     if numlabs>1
%         filename = ['Workspace1_Lab' num2str(labindex)];
%         save(filename);
%     else
%         filename = ['Workspace1_Lab' num2str(myid)];
%         load(filename);
%     end
% end


% (6) Create communication data for new MG level
% ==============================================

% bars_on_faces      = el2bar(:,13:24);
% new_nodes_on_faces = bars_on_faces + nVnod;
% nodface2sd             = zeros(nnod,1);
% for iel=1:nel_c
%     faces_c = Mesh_c.el2face(iel,:);
%     SD      = Mesh_c.face2sd(faces_c);
%     for iface=1:4
%         if SD(iface)
%             ind          = iface*3-2 : iface*3;
%             inod         = new_nodes_on_faces(iel,ind);
%             scatter3(GCOORD_SD(1,inod)',GCOORD_SD(2,inod)',GCOORD_SD(3,inod)',...
%                      700,'s','Linewidth',1,'MarkerEdgeColor','w');
%             nodface2sd(inod) = SD(iface);
%         end
%     end
% end
% nodface2sd1 = nodface2sd;

% Same as above but faster
nodface2sd = zeros(nnod,1);
for iface=1:4
    gface = el2face_c(:,iface);
    SD    = face2sd_c(gface);
    ind   = SD>0;
    nods  = el2bar(ind,10+iface*3:12+iface*3)+nVnod;
%     scatter3(GCOORD_SD(1,nods)',GCOORD_SD(2,nods)',GCOORD_SD(3,nods)',...
%              500,'o','Linewidth',2,'MarkerEdgeColor','y');
    nodface2sd(nods) = repmat(SD(ind),1,3);
end

% calculate pointer: element --> six edges
el2edge   = EL2NOD_SD(:,5:10)-nVnod;
EL2NOD_SD = EL2NOD_SD';
% calculate pointer: face --> three edges
face2edge = reshape(el2edge(:,[2 3 6 3 4 5 1 4 6 1 2 5])',3,[])';
[face2edge,~,ib] = unique_keep_order(face2edge);
nface     = size(face2edge,1);

% calculate pointer: element --> four faces
el2face = reshape(ib,4,nel)'; % element to face connectivity after doubles 
                              % have been merged (removed)

% addpath('Tools');
% patch_tetrahedron(60,GCOORD_SDc,EL2NOD_SDc,1,4,'k','-',1,14);
% axes_limits = [get(gca,'XLim') get(gca,'YLim') get(gca,'ZLim')];
% for iel=1:8
%     sfigure(60+iel);clf;
%     patch_tetrahedron(60+iel,GCOORD_SDc,EL2NOD_SDc,1,0,'k','-',1,14);
%     patch_tetrahedron(60+iel,GCOORD_SD,EL2NOD_SD,iel,4,'r','-',1,14);
%     axis(axes_limits);
% end

% Each parent element's face is split into 4 thus becomes the face of
% four child elements. The configuration is as follows:
% parent face, child element, child face
% ---------
% 1 , 2 , 1
% 1 , 3 , 1
% 1 , 4 , 1
% 1 , 6 , 4
% ---------
% 2 , 1 , 2
% 2 , 3 , 2
% 2 , 4 , 2
% 2 , 8 , 2
% ---------
% 3 , 1 , 3
% 3 , 2 , 3
% 3 , 4 , 3
% 3 , 5 , 3
% ---------
% 4 , 1 , 4
% 4 , 2 , 4
% 4 , 3 , 4
% 4 , 7 , 3

% Written in terms of the child elements. The number is the face of the
% child, the position (1..4) is the # of the face of the parent
% child 1 : [0 2 3 4]
% child 2 : [1 0 3 4]
% child 3 : [1 2 0 4]
% child 4 : [1 2 3 0]
% child 5 : [0 0 3 0] , i.e. face 3 of child 5 is on face 3 of the parent
% child 6 : [4 0 0 0] , i.e. face 4 of child 6 is on face 1 of the parent
% child 7 : [0 0 0 3] , i.e. face 3 of child 7 is on face 4 of the parent
% child 8 : [0 2 0 0] , i.e. face 2 of child 8 is on face 2 of the parent
%           position = face of parent (each column has 4 entires because
%                      each parent face is split into 4 triangles)
%           number = face of child
face2sd = zeros(nface,1);

% check all 1st childs and their faces 2,3, and 4
ind_els  = 1:8:nel;
faces_ch = [2 3 4]; % see list above
faces_pa = [2 3 4]; % see list above
for i=1:length(faces_ch)
    gface          = el2face(ind_els,faces_ch(i));
    gface_c        = el2face_c(:,faces_pa(i));
    face2sd(gface) = face2sd_c(gface_c);
end

% check all 2nd childs and their faces 1,3, and 4
ind_els  = 2:8:nel;
faces_ch = [1 3 4]; % see list above
faces_pa = [1 3 4]; % see list above
for i=1:length(faces_ch)
    gface          = el2face(ind_els,faces_ch(i));
    gface_c        = el2face_c(:,faces_pa(i));
    face2sd(gface) = face2sd_c(gface_c);
end

% check all 3rd childs and their faces 1,2, and 4
ind_els  = 3:8:nel;
faces_ch = [1 2 4]; % see list above
faces_pa = [1 2 4]; % see list above
for i=1:length(faces_ch)
    gface          = el2face(ind_els,faces_ch(i));
    gface_c        = el2face_c(:,faces_pa(i));
    face2sd(gface) = face2sd_c(gface_c);
end

% check all 4th childs and their faces 1,2, and 3
ind_els  = 4:8:nel;
faces_ch = [1 2 3]; % see list above
faces_pa = [1 2 3]; % see list above
for i=1:length(faces_ch)
    gface          = el2face(ind_els,faces_ch(i));
    gface_c        = el2face_c(:,faces_pa(i));
    face2sd(gface) = face2sd_c(gface_c);
end

% check all 5th childs and their face 3
ind_els        = 5:8:nel;
faces_ch       = 3; % see list above
faces_pa       = 3; % see list above
gface          = el2face(ind_els,faces_ch);
gface_c        = el2face_c(:,faces_pa);
face2sd(gface) = face2sd_c(gface_c);

% check all 6th childs and their face 4
ind_els        = 6:8:nel;
faces_ch       = 4; % see list above
faces_pa       = 1; % see list above
gface          = el2face(ind_els,faces_ch);
gface_c        = el2face_c(:,faces_pa);
face2sd(gface) = face2sd_c(gface_c);

% check all 7th childs and their face 3
ind_els        = 7:8:nel;
faces_ch       = 3; % see list above
faces_pa       = 4; % see list above
gface          = el2face(ind_els,faces_ch);
gface_c        = el2face_c(:,faces_pa);
face2sd(gface) = face2sd_c(gface_c);

% check all 8th childs and their face 2
ind_els        = 8:8:nel;
faces_ch       = 2; % see list above
faces_pa       = 2; % see list above
gface          = el2face(ind_els,faces_ch);
gface_c        = el2face_c(:,faces_pa);
face2sd(gface) = face2sd_c(gface_c);


% Edges that are on edges of parent elements
% 1st child
els    = 1:8:nel;   %  5  8  9
edges1 = el2edge(els,[ 1  4  5]);
% 2nd child
els    = 2:8:nel;   %  5  6 10
edges2 = el2edge(els,[ 1  2  6]);
% 3rd child
els    = 3:8:nel;   %  6  7  9
edges3 = el2edge(els,[ 2  3  5]);
% 4th child
els    = 4:8:nel;   %  7  8 10
edges4 = el2edge(els,[ 3  4  6]);

nedge  = nnod-nVnod;
edgesA = unique([edges1 edges2 edges3 edges4]);
edgesB = setdiff(1:nedge,edgesA);
if length(edgesA)~=2*length(edge2sd_c)
    error('Edge calculation failed');
end

% Generate new pointer edge to subdomain
edge2sd  = zeros(nedge,size(edge2sd_c,2));
nNB_edge = zeros(nedge,1);

% Initialize a unique list of nodes (i.e. if all SDs merge there unique
% nodes, none is missing and none is counted twice; required for parallel
% inner product)
unique_nodes = 1:nnod;
mynod_SDBc   = COMM.mynod_SDB{img};
for iNB=1:COMM.nNB
    labBarrier
    NB = COMM.NB(iNB);
    if ~NB
        continue
    end
    
    % Find fine mesh vertex nodes (i.e. coarse mesh nodes) that are in the
    % list of shared nodes (with neighbor NB)
    mynods       = mynod_SDBc{iNB};
    [tf,~]       = ismember(pointer_Edge2Vnod,mynods); % find mynods in pointer_Edge2Vnod
    tf(edgesB,:) = 0;
    myEdges      = find(sum(tf,2)==2)';

    % Create a pointer element edge to neighboring subdomain
    nNB_edge(myEdges) = nNB_edge(myEdges) + 1;
    ind          = sub2ind([nedge nsd],myEdges(:),nNB_edge(myEdges));
    edge2sd(ind) = NB;                  % write edge to SD pointer
    
    myEdgeNods1 = nVnod + myEdges;
    myEdgeNods2 = find(nodface2sd==NB);
    
    if FigNo
        sfigure(FigNo);
        sz = 30 + (iNB-1)*100;
        sc = ['g','m','b'];
        scatter3(GCOORD_SD(1,mynods)',GCOORD_SD(2,mynods)',GCOORD_SD(3,mynods)',...
                 sz,'o','Linewidth',1.5,'MarkerEdgeColor',sc(iNB));
        scatter3(GCOORD_SD(1,myEdgeNods1)',GCOORD_SD(2,myEdgeNods1)',GCOORD_SD(3,myEdgeNods1)',...
                 sz,'o','Linewidth',1.5,'MarkerEdgeColor',sc(iNB));
        scatter3(GCOORD_SD(1,myEdgeNods2)',GCOORD_SD(2,myEdgeNods2)',GCOORD_SD(3,myEdgeNods2)',...
                 sz,'o','Linewidth',2,'MarkerEdgeColor',sc(iNB));
        if iNB==COMM.nNB
            saveas(gcf,sprintf('./Mesh_MG%1i_SD%1ix%1i',img+1,myid,nsd),'fig');
        end
    end
    
    myEdgeNods = union(myEdgeNods1,myEdgeNods2);
    
    if numlabs>1
        % Bring new edge nodes into the right order so they match the
        % neighbor's nodes. Use node coordinates to do so.
        my_xyz = single(GCOORD_SD(:,myEdgeNods))';
        nb_xyz = COMM.sendnrecv_1NB(NB,my_xyz); % *COMMUNICATION*
    
        if myid<NB % Only one SD must re-order!
            [tf,perm] = ismember(nb_xyz,my_xyz,'rows');
            if any(~tf)
                error('Generation of subdomain boundary failed.');
            end
            % re-order coordinates and check that the match
            my_xyz = my_xyz(perm,:);
            if max(max(abs(my_xyz-nb_xyz)))>1e-8
                error('Generation of subdomain boundary failed.');
            end
            % re-order new edge nodes
            myEdgeNods = myEdgeNods(perm);
        end
    end
    
%     if img==2
%         if numlabs>1
%             filename = ['Workspace2_Lab' num2str(labindex)];
%             save(filename);
%             stop
%         else
%             my_xyz=[]; nb_xyz=[];
%             filename = ['Workspace2_Lab' num2str(myid)];
%             load(filename);
%         end
%     end
    
    % all nodes shared with SD "NB"
    mynod_SDB                  = [mynods myEdgeNods];
    COMM.nnod_SDB {img-1}(iNB) = length(mynod_SDB);
    COMM.mynod_SDB{img-1}{iNB} = mynod_SDB;
    
    % Pressure nodes in SD boundary
    myPdof_SDB                  = mynod_SDB(mynod_SDB<=nVnod);
    COMM.nPdof_SDB{img-1}(iNB)  = length(myPdof_SDB);
    COMM.myPdof_SDB{img-1}{iNB} = myPdof_SDB;

    % Velocity dofs on SD boundary
    myUdof_SDB                  = nod2dof(mynod_SDB,3);
%     myUdof_SDB                = [3*mynod_SDB(:)'-2;
%                                  3*mynod_SDB(:)'-1;
%                                  3*mynod_SDB(:)'   ];
    myUdof_SDB                  = myUdof_SDB(:)';
    COMM.nUdof_SDB{img-1}(iNB)  = length(myUdof_SDB);
    COMM.myUdof_SDB{img-1}{iNB} = myUdof_SDB;
    
    % Update unique node list (nodes on a SD boundary are assigned to the
    % SD with higher index)
    if myid<NB
        unique_nodes(mynod_SDB) = 0;
    end
end

% Write unique lists for temperature and pressure nodes, and velocity dofs
unique_nodes             = uint32(find(unique_nodes));
COMM.unique_nodes{img-1} = unique_nodes;
COMM.unique_Udofs{img-1} = nod2dof(unique_nodes,3);
COMM.unique_Pdofs{img-1} = unique_nodes( unique_nodes<=nVnod );

end % END OF SUBFUNCTION generate_next_MG_level_old

% #########################################################################

function COMM = communication_function_handles(COMM)

% Handles of all communication functions used in the code
COMM.sendnrecv_1NB = @sendnrecv_1NB; % uses labSendReceive
COMM.send_1NB      = @send_1NB;
COMM.recv_1NB      = @recv_1NB;
% COMM.sum_all       = @sum_all;
COMM.sum_all       = @sum_all_v2; % alternative way of doing the parallel sum
COMM.dot           = @dot_p;
% COMM.dot           = @dot_p_v2; % alternative way of doing the parallel dot
COMM.norm2         = @norm2_p;
COMM.normdf        = @normdf_p;   % alternative norm
COMM.minLabs       = @minLabs;
COMM.maxLabs       = @maxLabs;
COMM.vcat          = @vcat_p;
COMM.hcat          = @hcat_p;
COMM.unique        = @unique_p;
COMM.sum_SDB       = @sum_SDB; % uses labSendReceive

end % END OF SUBFUNCTION communication_function_handles

% #########################################################################

function NBdata = sendnrecv_1NB(NB,mydata,tag)
    if nargin==2
        tag=1;
    end
    if NB<0 || NB>numlabs
        NBdata = nan(size(mydata)); % to make it run in serial mode for debugging
        return
    else
        NBdata = labSendReceive( NB , NB , mydata , tag );
    end
end % END OF SUBFUNCTION sendnrecv_1NB

% #########################################################################

function send_1NB(NB,mydata,tag)
    if nargin==2
        tag=1;
    end
    if NB<0 || NB>numlabs
        return
    else
    %     disp(['Lab ' num2str(labindex) ' sends to NB ' num2str(NB)]);
        labSend( mydata , NB , tag );
    end
end % END OF SUBFUNCTION send_1NB

% #########################################################################

function NBdata = recv_1NB(NB,tag)
    if nargin==1
        tag=1;
    end
    if NB<0 || NB>numlabs
        NBdata = nan; % to make it run in serial mode for debugging
        return
    else
    %     disp(['Lab ' num2str(labindex) ' receives from NB ' num2str(NB)]);
        NBdata = labReceive( NB , tag );
    end
end % END OF SUBFUNCTION recv_1NB

% #########################################################################

function sum = sum_all(mydata)
    % disp('summing up');
    sum = 0;                                     
    for ilab=1:numlabs
        labBarrier
        if ilab==labindex
            sum = sum + labBroadcast( labindex , mydata );
        else
            sum = sum + labBroadcast( ilab );
        end
    end
    % disp('done');
end % END OF SUBFUNCTION sum_all

% #########################################################################

function sum = sum_all_v2(mydata)
    sum = mydata;
    % lab 1 sends to 2, 2 adds it's part, 2 sends to 3, ...
    for ilab=1:numlabs-1
        if labindex==ilab
            labSend( sum , ilab+1 );
        elseif labindex==ilab+1
            sum = sum + labReceive( ilab );
        end 
    end
    % last lab has the complete dot product and sends it to all others
    if labindex==numlabs
        labBroadcast( labindex , sum );
    else
        sum = labBroadcast( numlabs );
    end
end % END OF SUBFUNCTION sum_all_v2

% #########################################################################

function ab = dot_p(a,b,unique_dofs)
    % disp('dot product');
    ab_all           = zeros(numlabs,1);
    ab_all(labindex) = dot( a(unique_dofs) , b(unique_dofs) );
    for ilab=1:numlabs
        labBarrier
        if labindex==ilab
            ab_all(ilab) = labBroadcast( ilab , ab_all(ilab) );
        else
            ab_all(ilab) = labBroadcast( ilab );
        end
    end
    ab = sum(ab_all); % sum up 'ab' over the domain
    % disp('done');
end % END OF SUBFUNCTION dot_p

% #########################################################################

function ab = dot_p_v2(a,b,unique_dofs)
    % each lab does it's part of the dot product
    ab = dot( a(unique_dofs) , b(unique_dofs) );
    
    % lab 1 sends to 2, 2 adds it's part, 2 sends to 3, ...
    for ilab=1:numlabs-1
        if labindex==ilab
            labSend( ab , ilab+1 );
        elseif labindex==ilab+1
            ab = ab + labReceive( ilab );
        end 
    end
    
    % last lab has the complete dot product and sends it to all others
    if labindex==numlabs
        labBroadcast( labindex , ab );
    else
        ab = labBroadcast( numlabs );
    end
end % END OF SUBFUNCTION dot_p_v2

% #########################################################################

function rms = norm2_p(a,unique_dofs)
    rms = dot_p(a,a,unique_dofs);
    rms = sqrt(rms);
end % END OF SUBFUNCTION norm2_p

% #########################################################################

function rms = normdf_p(a,unique_dofs)
    rms = dot_p(a,a,unique_dofs);
    nD  = sum_all(length(unique_dofs));
    rms = sqrt(rms/nD);
end % END OF SUBFUNCTION normdf_p

% #########################################################################

function value = maxLabs(myvalue)
    myvalue = max(myvalue);
    value   = myvalue;
    for ilab=1:numlabs
        labBarrier
        if labindex==ilab
            labBroadcast( ilab , myvalue );
        else
            value = max(value,labBroadcast( ilab ));
        end
    end
end % END OF SUBFUNCTION maxLabs

% #########################################################################

function value = minLabs(myvalue)
    myvalue = min(myvalue);
    value   = myvalue;
    for ilab=1:numlabs
        labBarrier
        if labindex==ilab
            labBroadcast( ilab , myvalue );
        else
            value = min(value,labBroadcast( ilab ));
        end
    end
end % END OF SUBFUNCTION minLabs

% #########################################################################

function a = vcat_p(mydata)
    a = mydata;
    for ilab=1:numlabs
        labBarrier
        if ilab==labindex
            labBroadcast( labindex , mydata );
        else
            a = [a; labBroadcast( ilab )]; %#ok<AGROW>
        end
    end
end % END OF SUBFUNCTION vcat_p

% #########################################################################

function a = hcat_p(mydata)
    a = mydata;
    for ilab=1:numlabs
        labBarrier
        if ilab==labindex
            labBroadcast( labindex , mydata );
        else
            a = [a labBroadcast( ilab )]; %#ok<AGROW>
        end
    end
end % END OF SUBFUNCTION hcat_p

% #########################################################################

function a = unique_p(mydata)
    a = mydata(:);
    for ilab=1:numlabs
        labBarrier
        if ilab==labindex
            labBroadcast( labindex , mydata );
        else
            b = labBroadcast( ilab );
            a = [a; b(:)]; %#ok<AGROW>
        end
    end
    a = unique(a);
end % END OF SUBFUNCTION unique_p

% #########################################################################

function data = sum_SDB(mydata,NBs,mydofs_allNB)
    % disp('sumU_SD_Bound2');
    data = mydata; % otherwise summed up data is sent to neighbors 
                   % causing multiple summations
    for iNB=1:length(NBs)
%         labBarrier
        if NBs(iNB)~=0
            mydofs = mydofs_allNB{iNB};
            NBdata = labSendReceive( NBs(iNB) , NBs(iNB) , mydata(mydofs) );
            data(mydofs) = data(mydofs) + NBdata;
        end
    end
    % disp('done');
end % END OF SUBFUNCTION sum_SDB