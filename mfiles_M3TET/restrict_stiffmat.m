function [KK,Ic2f_U] = restrict_stiffmat(MESH,K1,ifree,fid_log)
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

fprintf(fid_log,' Restricting stiffness matrix to all MG levels...');

nmg           = MESH.nmg; % number of geometric multigrid levels
nmesh         = length(MESH.EL2NOD); % total number of meshes (there can be an additional, non-nested base mesh)
nUdof_mg      = zeros(nmesh,1);
nUdof_mg_free = zeros(nmesh,1);
for img=1:nmesh
    nUdof_mg(img) = 3*max(max(MESH.EL2NOD{img}));
    if img>nmg
        nUdof_mg_free(img) = nUdof_mg(img);
    else
        nUdof_mg_free(img) = length(find(ifree<=nUdof_mg(img)));
    end
end

ndofblo = 100000; % number of rows/columns that are multiplied at once
KK      = cell(nmesh,1);   % allocate cell-array for stiffness matrices
KK{1}   = full(diag(K1));  % store the diagonal of K1 in the cell-array KK
                           % matrix K1 is kept separately for performance
Ic2f_U  = cell(nmesh-1,1); % allocate cell-array for interpolation matrices

% Recursively write interpolation and stiffness matrices into cell arrays
for img=1:nmg-1
    iUfree_c    = intersect(1:nUdof_mg(img+1),ifree);
    Ic2f_U{img} = MESH.Ic2f{img}(ifree,iUfree_c);
    ifree       = iUfree_c;
    if img==1
        KK{img+1} = restrict_matrix_blockwise(Ic2f_U{img},K1,ndofblo);
    else
        KK{img+1} = Ic2f_U{img}' * KK{img} * Ic2f_U{img};
    end
end

if nmesh>nmg
    Ic2f_U{nmg} = MESH.Ic2f{nmg}(ifree,:);
    KK{nmg+1}   = Ic2f_U{nmg}' * KK{nmg} * Ic2f_U{nmg};
end

fprintf(fid_log,'done.\n');
fprintf(fid_log,' Problem size (nUdof)                        = ');
fprintf(fid_log,'%1i  ',nUdof_mg);
fprintf(fid_log,'\n');
fprintf(fid_log,' Number of free Udofs                        = ');
fprintf(fid_log,'%1i  ',nUdof_mg_free);
fprintf(fid_log,'\n');

end % END OF SUBFUNCTION restrict_stiffmat_freedofs

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function K2 = restrict_matrix_blockwise(Ic2f,K1,ndofblo)

% Purpose: Block-wise multiplication of sparse matrices to reduce the peak
%          in memory usage that occurs when MATLAB re-defines the sparse
%          structure of the resulting matrix.
%          Here the restriction K2 = R * K1 * I is performed, where R=I' is
%          the restriction operator and I the interpolation operator
% Author:  Joerg Hasenclever, IFG, University of Hamburg
%          now at: GEOMAR Helmholtz Centre for Ocean Science, Kiel
ndof    = length(K1);
ndofblo = min(ndof,ndofblo);
nblo    = ceil(ndof/ndofblo);
ia      = 1;
ib      = ndofblo;
Restrc  = Ic2f';
if ~any(any(full(triu(K1(1:3,1:3),1))))
    KD1 = full(diag(K1));
    K1  = tril(K1,-1);
    for ii=1:nblo
        % Lower triangle
        tmp1 = Restrc * K1(:,ia:ib);
        if ii==1
            K2   = tmp1*Ic2f(ia:ib,:);
        else
            tmp2 = tmp1*Ic2f(ia:ib,:);
            K2   = K2 + tmp2;
        end    

        % Upper triangle
        tmp1 = Restrc(:,1:ib)*K1(ia:ib,1:ib)';
        tmp2 = tmp1*Ic2f(ia:ib,:);
        K2   = K2 + tmp2;

        % ----------------
        ia  = ib+1;
        ib  = min(ib+ndofblo,ndof);
    end
    % Diagonal
    tmp2 = Restrc * diag(sparse(KD1)) * Ic2f;
    K2   = K2 + tmp2;
    
else
    for ii=1:nblo
        if ii==1
            K2   = Restrc*K1(:,ia:ib)*Ic2f(ia:ib,:);
        else
            tmp2 = Restrc*K1(:,ia:ib)*Ic2f(ia:ib,:);
            K2   = K2 + tmp2;
        end
        ia = ib+1;
        ib = min(ib+ndofblo,ndof);
    end
end

end % END OF SUBFUNCTION restrict_matrix_blockwise