function addpaths_suitesparse()

% Check if SuiteSparse functions are already accessible:
if exist('sparse2','file') && exist('cs_lsolve','file') && exist('cs_chol','file')
    return % functions have been located
end

% Nope, we have to add paths for this computer:
hostname = get_hostname;
switch hostname
    case 'b3pc14'
        return % path already known by Matlab
    case 'b3pc16'
        return % path already known by Matlab
    case 'Maupiti2'
        return % path already known by Matlab
    case 'snaefell'
        path2suitesparse = '/Users/joha/Work/Numerical_Libraries/SuiteSparse/';
    case 'fuego'
        path2suitesparse = '/scratch/local1/spielwiese/SuiteSparse';
    case 'santiaguito'
        path2suitesparse = '/scratch/local1/spielwiese/SuiteSparse';
    otherwise
        if strcmp(hostname(1:4),'node') || strcmp(hostname,'master')
            path2suitesparse = '/home/jhasenclever/Matlab_Tools/SuiteSparse_unix64/';
        elseif strcmp(hostname(1:4),'rzcl')
            path2suitesparse = '/work_j/smomw219/Matlab_Tools/SuiteSparse_unix64/';
        else
            error(' SuiteSparse location not defined for this computer.\n You have to edit file "addpaths_suitesparse.m"!');
        end
end

addpath([path2suitesparse '']);
addpath([path2suitesparse '/UMFPACK/MATLAB']);
addpath([path2suitesparse '/CHOLMOD/MATLAB']);
addpath([path2suitesparse '/AMD/MATLAB']);
addpath([path2suitesparse '/COLAMD/MATLAB']);
addpath([path2suitesparse '/CCOLAMD/MATLAB']);
addpath([path2suitesparse '/CAMD/MATLAB']);
addpath([path2suitesparse '/CXSparse/MATLAB/UFget']);
addpath([path2suitesparse '/CXSparse/MATLAB/Demo']);
addpath([path2suitesparse '/CXSparse/MATLAB/CSparse']);
addpath([path2suitesparse '/LDL/MATLAB']);
addpath([path2suitesparse '/BTF/MATLAB']);
addpath([path2suitesparse '/KLU/MATLAB']);
addpath([path2suitesparse '/SPQR/MATLAB']);
addpath([path2suitesparse '/RBio/RBio']);
addpath([path2suitesparse '/MATLAB_Tools']);
addpath([path2suitesparse '/MATLAB_Tools/Factorize']);
addpath([path2suitesparse '/MATLAB_Tools/MESHND']);
addpath([path2suitesparse '/MATLAB_Tools/LINFACTOR']);
addpath([path2suitesparse '/MATLAB_Tools/find_components']);
addpath([path2suitesparse '/MATLAB_Tools/GEE']);
addpath([path2suitesparse '/MATLAB_Tools/shellgui']);
addpath([path2suitesparse '/MATLAB_Tools/waitmex']);
addpath([path2suitesparse '/MATLAB_Tools/spqr_rank']);
addpath([path2suitesparse '/MATLAB_Tools/spqr_rank/SJget']);
addpath([path2suitesparse '/MATLAB_Tools/UFcollection']);
addpath([path2suitesparse '/MATLAB_Tools/SSMULT']);
addpath([path2suitesparse '/MATLAB_Tools/dimacs10']);
addpath([path2suitesparse '/MATLAB_Tools/spok']);
addpath([path2suitesparse '/MATLAB_Tools/sparseinv']);


% Check that suitesparse functions are accessible now
try
    A = eye(5); test = 1;
    A = sparse2(A); test = 2;
    A = cs_transpose(A); test = 3;
    L = cs_chol(A); test = 4;
    x = cs_ltsolve(L,cs_lsolve(L,ones(5,1)));
catch
    fid = fopen(['Lab' num2str(numlabs,2) 'x' num2str_d(labindex,2) '_SuiteSparse_ERROR.log'],'w');
    fprintf(fid,'Machine "%s"\n',hostname);
    fprintf(fid,'Paths to SuiteSparse: "%s"\n',path2suitesparse);
    switch test
        case 1
            err = 'sparse2';
        case 2
            err = 'cs_transpose';
        case 3
            err = 'cs_chol';
        case 4
            err = 'cs_ltsolve and/or cs_lsolve';
    end
    fprintf(fid,'\n ERROR: Failed to run command "%s"\n',err);
    p = path; fprintf(fid,'All paths:\n %s',p);
    fclose(fid);
    error(' Failed to run command "%s" after adding paths to SuiteSparse!',err);
end

end % END OF FUNCTION addpaths_suitesparse