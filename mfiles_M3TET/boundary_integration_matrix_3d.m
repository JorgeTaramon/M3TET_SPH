function NN_INTEG = boundary_integration_matrix_3d(GCOORD,EL2NOD,nods_bnd)

% Calculates the matrix NN_INTEG that, when multiplied a nodal vector, will
%  integrate these nodal values along the domain boundary
switch size(EL2NOD,1)
    case 4
        nnod_face = 3; % nodes on element segment
        nip       = 3;
    case 10
        nnod_face = 6; % nodes on element segment
        nip       = 3;
end

NN_INTEG = assemble_integration_matrix_tetra(EL2NOD,GCOORD,nods_bnd,nnod_face,nip);

end % END OF FUNCTION boundary_integration_matrix_3d

% #########################################################################
%                               SUBFUNCTIONS
% #########################################################################

function NN_INTEG = assemble_integration_matrix_tetra(EL2NOD,GCOORD,nods_bnd,nnod_face,nip)

switch nip
    case 1
        % ONE POINT INTEGRATION
        r  = 1/3;
        s  = 1/3;
        wt = 1;
    case 3
        % THREE POINT INTEGRATION (version 1)
        % Taken from Zienkiewicz Vol 1 (4th) edition p.176 (also in Hughes p.173)
        r  = [0.5 0.5 0.0]; % r
        s  = [0.5 0.0 0.5]; % s
        wt = [1/3 1/3 1/3]; % weights for integration (sum up to one)
        
%         % THREE POINT INTEGRATION (version 2)
%         r  = [1/6 2/3 1/6]; % r
%         s  = [1/6 1/6 2/3]; % s
%         wt = [1/3 1/3 1/3]; % weights for integration (sum up to one)
    case 7
        % SEVEN POINT INTEGRATION
        % Taken from Zienkiewicz Vol 1 (4th) edition p.177 (also in Hughes p.173)
        % Integrates 5th-order polynomials exactly.
        a  = 0.470142064105115; b = 0.059715871789770;
        c  = 0.101286507323456; d = 0.797426985353087;
        w1 = 0.225; w2 = 0.132394152788506; w3 = 0.125939180544827;
        r  = [1/3 a b a c d c]; % r of all points
        s  = [1/3 a a b c c d]; % s of all points
        wt = [w1 w2 w2 w2 w3 w3 w3]; % integration weights at each point

    otherwise
        error('nip must be 1,3, or 7');
end
x_ip   = [r;s];
[N,dN] = sf_dsf_tri367(x_ip,nnod_face,'cell');

% TEMPERATURE BOTTOM FLUX BC
bc_els = find(sum(ismember(EL2NOD,nods_bnd),1)==nnod_face); % find elements with nnod_face points on boundary

% storage
kki = zeros(nnod_face*nnod_face,length(bc_els)); % storage for big sparse matrix
kkj = zeros(nnod_face*nnod_face,length(bc_els));
kk1 = zeros(nnod_face*nnod_face,length(bc_els));

% Assemble the matrix that multiplied onto the prescribed flux BC's
% gives the integrated flux values that get added to the rhs of the
% temperature equation
for ii=1:length(bc_els)
    iel        = bc_els(ii); % global element number
    EL2NOD_tri = tetra_face_connectivity(GCOORD,EL2NOD(:,iel),nods_bnd); % *SUBFUNCTION*
    
    % Calculate area of tetrahedron face
    xyz_el     = GCOORD(:,EL2NOD_tri);
    a          = sqrt(sum( (xyz_el(:,1)-xyz_el(:,2)).^2 ));
    b          = sqrt(sum( (xyz_el(:,2)-xyz_el(:,3)).^2 ));
    c          = sqrt(sum( (xyz_el(:,3)-xyz_el(:,1)).^2 ));
    d          = 0.5*(a+b+c); % semiperimeter (half of triangle's perimeter)
    area       = sqrt( d*(d-a)*(d-b)*(d-c) ); % Heron's formula
    
    if abs(min(xyz_el(3,:))-max(xyz_el(3,:))) < 1e-8
        J         = xyz_el(1:2,:) * dN{1}';
        if (2*area - abs(det(J))) > 1e-12*area
            error('Area calculated using Herons formula should be equal to 2*det(J).');
        end
    end
    
    f = zeros(nnod_face);            % local element matrix
    for ip=1:length(wt)
        f = f + (N{ip}*N{ip}') * area * wt(ip);
    end

    rows      = double(EL2NOD_tri(:))*ones(1,nnod_face);
    cols      = ones(nnod_face,1)*double(EL2NOD_tri(:))';
    kk1(:,ii) = f(:);
    kki(:,ii) = rows(:);
    kkj(:,ii) = cols(:);
end

nnod     = max(EL2NOD(:));
NN_INTEG = sparse3(kki(:),kkj(:),kk1(:),nnod,nnod); % create big sparse matrix

end % END OF SUBFUNCTION assemble_integration_matrix_tetra

% #########################################################################

function EL2NOD_tri = tetra_face_connectivity(GCOORD,EL2NOD_tet,nods_bnd)

[EL2NOD_tri,~,~] = intersect(EL2NOD_tet,nods_bnd);
nnod_face        = length(EL2NOD_tri);

% Make sure numbering is counterclockwise
a   = GCOORD(:,EL2NOD_tri(2)) - GCOORD(:,EL2NOD_tri(1));
b   = GCOORD(:,EL2NOD_tri(3)) - GCOORD(:,EL2NOD_tri(1));
axb = cross(a,b);
if axb(3)<0
    EL2NOD_tri(1:3) = EL2NOD_tri([1 3 2]);
end

if nnod_face==6
    xyz_el = GCOORD(:,EL2NOD_tri);
%     figure(20);clf;
%     for i=1:nnod_face
%         scatter(xyz_el(1,i),xyz_el(2,i),500); hold on
%         text(xyz_el(1,i),xyz_el(2,i),num2str(i),'Fontsize',18,'Fontweight','bold');
%     end
    
    xyz4   = mean(xyz_el(:,[1 2]),2);
    xyz5   = mean(xyz_el(:,[2 3]),2);
    xyz6   = mean(xyz_el(:,[1 3]),2);
    ind4   = 3 + find( sum(abs(xyz_el(:,4:6)-repmat(xyz4,1,3))) < 1e-8 );
    ind5   = 3 + find( sum(abs(xyz_el(:,4:6)-repmat(xyz5,1,3))) < 1e-8 );
    ind6   = 3 + find( sum(abs(xyz_el(:,4:6)-repmat(xyz6,1,3))) < 1e-8 );
    EL2NOD_tri(4:6) = EL2NOD_tri([ind4 ind5 ind6]);
    
%     xyz_el = GCOORD(:,EL2NOD_tri);
%     figure(20);clf;
%     for i=1:nnod_face
%         scatter(xyz_el(1,i),xyz_el(2,i),500); hold on
%         text(xyz_el(1,i),xyz_el(2,i),num2str(i),'Fontsize',18,'Fontweight','bold');
%     end
end

end % END OF SUBFUNCTION tetra_face_connectivity

% #########################################################################

function NN_INTEG = assemble_integration_matrix_hexa(EL2NOD,GCOORD,nods_bnd,nnod_face,nip)

nip_1D = sqrt(nip);
switch nip_1D
    case 1
        pt_1D = 0;
        wt_1D = 2;
        
    case 2
        pt_1D = [-1/sqrt(3) 1/sqrt(3)];
        wt_1D = [1 1];
        
    case 3
        pt_1D = [-sqrt(3/5) 0 sqrt(3/5)];
        wt_1D = [5/9 8/9 5/9];
        
    otherwise
        error(' nip must be 1,4,9,...');
end
r  = zeros(1,nip);
s  = zeros(1,nip);
wt = zeros(1,nip);
i  = 0;
for i1=1:nip_1D
    for i2=1:nip_1D
        i     = i+1;
        r(i)  = pt_1D(i1);
        s(i)  = pt_1D(i2);
        wt(i) = wt_1D(i1) * wt_1D(i2);
    end
end
wt = 1 * wt ./ sum(wt);

nnod   = max(EL2NOD(:));
[N,dN] = shapefct_quad(nnod_face,[r(:)';,s(:)']);

% TEMPERATURE BOTTOM FLUX BC
bc_els = find(sum(ismember(EL2NOD,nods_bnd),1)==nnod_face); % find elements with nnod_face points on boundary

% storage
kki = zeros(nnod_face*nnod_face,length(bc_els)); % storage for big sparse matrix
kkj = zeros(nnod_face*nnod_face,length(bc_els));
kk1 = zeros(nnod_face*nnod_face,length(bc_els));

% Assemble the matrix that multiplied onto the prescribed flux BC's
% gives the integrated flux values that get added to the rhs of the
% temperature equation
for ii=1:length(bc_els)
    iel      = bc_els(ii); % global element number    
    [~,p1,~] = intersect(EL2NOD(:,iel),nods_bnd);
    
    % Need to resort the nodes on the element's face so that they are in
    % counter-clockwise order when viewed from outside the domain.
    if     isempty(setdiff(p1,[1 2 3 4])) % FACE #1
        p1 = [4 3 2 1];
    elseif isempty(setdiff(p1,[1 2 5 6])) % FACE #2
        p1 = [1 2 6 5];
    elseif isempty(setdiff(p1,[2 3 6 7])) % FACE #3
        p1 = [2 3 7 6];
    elseif isempty(setdiff(p1,[3 4 7 8])) % FACE #4
        p1 = [3 4 8 7];
    elseif isempty(setdiff(p1,[1 4 5 8])) % FACE #5
        p1 = [4 1 5 8];
    elseif isempty(setdiff(p1,[5 6 7 8])) % FACE #6
        p1 = [5 6 7 8];
    end
    p1      = p1(end:-1:1); % WHY ??????
    nod_glb = EL2NOD(p1,iel);
    
    % Resort the 3 node on the tetra face so that they are in
    % counter-clockwise order (like in a 2D FE triangle code)
    if nnod_face==8
        error('2B coded');
    end

    % Calculate area of quadrilateral face
    xyz_el = GCOORD(:,nod_glb);
    J      = xyz_el(1:2,:) * dN(:,:,1)';
    detJ   = det(J);
    if detJ<0
        error('Negative Jacobian');
    end
    area   = 4*detJ;
    f      = zeros(nnod_face);            % local element matrix
    for ip=1:length(wt)
        f = f + (N(ip,:)'*N(ip,:)) * area * wt(ip);
    end

    rows      = double(nod_glb)*ones(1,nnod_face);
    cols      = ones(nnod_face,1)*double(nod_glb)';
    kk1(:,ii) = f(:);
    kki(:,ii) = rows(:);
    kkj(:,ii) = cols(:);
end

NN_INTEG = sparse3(kki(:),kkj(:),kk1(:),nnod,nnod); % create big sparse matrix

end % END OF SUBFUNCTION assemble_integration_matrix_hexa