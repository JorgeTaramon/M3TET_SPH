function addpaths_mutils()

% Check if MUTILS functions are already accessible:
if exist('tsearch2','file') && exist('metis','file') && ...
   exist('lchol','file') && exist('sparse_create','file')
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
        path2mutils      = '/Users/joha/Work/Numerical_Libraries/mutils-0.4-1/';
    case 'fuego'
        return % path already known by Matlab
    case 'colima'
        return % path already known by Matlab
    case 'jorgeMAC'
        path2mutils      = '/Users/jorge/Dropbox/GEOMAR/mutils-0.4-2';
    case 'clusterRHUL'
        path2mutils      = '/home/mat1/mutils-0.4-2_par';
    case 'jorgeLAPTOP'
        path2mutils      = 'C:\Users\Tara\Documents\MATLAB\mutils-0.4-2';
    case 'jmc1'
        path2mutils      = '/home/mat1/mutils-0.4-2_par';
    case 'clusterRHUL_2'
        path2mutils      = '/home/mat2/mutils-0.4-2_par';
    case 'glmsc12'
        path2mutils      = '/home/mat2/mutils-0.4-2_par';
    otherwise
        if strcmp(hostname(1:4),'node') || strcmp(hostname,'master')
            path2mutils      = '/home/jhasenclever/Matlab_Tools/mutils-0.4-2';
        elseif strcmp(hostname(1:4),'rzcl')
            path2mutils      = '/work_j/smomw219/Matlab_Tools/mutils-0.4-2/';
        else
            error(' MUTILS location not defined for this computer.\n You have to edit file "addpaths_mutils.m"!');
        end
end

if ~isempty(path2mutils)
    addpath( path2mutils );
    addpath([path2mutils '/triangle']);
    addpath([path2mutils '/SuiteSparse']);
    addpath([path2mutils '/mutils']);
    addpath([path2mutils '/mutils/quadtree']);
    addpath([path2mutils '/mutils/interp']);
    addpath([path2mutils '/mutils/reorder']);
    addpath([path2mutils '/mutils/sparse']);
else
    error(' Please install MUTILS on this computer and add the path in this function.');
end

end % END OF FUNCTION addpaths_mutils