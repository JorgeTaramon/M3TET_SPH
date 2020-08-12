function OPTS = DEFINE_OPTS(config,nlab,queue)

if nargin<2
    error('call "OPTS = DEFINE_OPTS(config,nlab)"');
end

if strcmp(config,'local')
    OPTS.config      = 'local';
    OPTS.nlab        = nlab;
    OPTS.path2source = pwd;
    OPTS.matlab_version = version('-release');
    if isunix
        if strcmp(get_hostname,'fuego') || strcmp(get_hostname,'colima')
            OPTS.path2mutils = '/opt/Matlab-Libs/mutils-0.4-2_with_openmp';
            OPTS.path2tools  = '/opt/Matlab-Libs/tools';
        elseif strcmp(get_hostname,'nesh-fe')
            OPTS.path2mutils = '~/MATLAB/Matlab-Libs/mutils-0.4-2_with_openmp';
            OPTS.path2tools  = '~/MATLAB/Matlab-Libs/tools';
        elseif strcmp(get_hostname,'jorgeMAC')
            OPTS.path2mutils = '/Users/jorge/Dropbox/GEOMAR/mutils-0.4-2';
            OPTS.path2tools  = [];
        elseif strcmp(get_hostname,'clusterRHUL')
            OPTS.path2mutils = '/home/mat1/mutils-0.4-2_par';
            OPTS.path2tools  = [];
        elseif strcmp(get_hostname,'jmc1')
            OPTS.path2mutils = '/home/mat1/mutils-0.4-2_par';
            OPTS.path2tools  = [];
        else
            disp(' You have to define "OPTS.path2mutils" and "OPTS.path2tools" manually.');
        end
    else
        OPTS.path2mutils = 'C:\Matlab-Libs\mutils-0.4-2_openmp_new';
        OPTS.path2tools  = 'C:\Matlab-Libs\tools';
    end
    
elseif strcmp(config,'nec')
    if nargin<3
        error('call "OPTS = DEFINE_OPTS(config,nlab,queue)"');
    end
    % Use this function to quickly define structure "OPTS" required
    % for a parallel job submission
    OPTS.config   = 'nec';

    % batch classes: clexpress, clmedium, cllong
    switch queue
        case 'clexpress'
            if nlab>2
                error('nlab<=2 for "clexpress"');
            end
            OPTS.queue    = queue;
            OPTS.walltime = '2:00:00';  % max: 2h
            OPTS.nlab     = nlab;       % total: 2

        case 'clmedium'
            OPTS.queue    = queue; 
            OPTS.walltime = '48:00:00'; % max: 48h
            OPTS.nlab     = nlab;       % total: 78

        case 'cllong'
            OPTS.queue    = queue; 
            OPTS.walltime = '100:00:00';% max: 100h
            OPTS.nlab     = nlab;       % total: 30

        case 'clfocean'
            if nlab>4
                error('nlab<=4 for "clfocean"');
            end
            OPTS.queue    = queue; 
            OPTS.walltime = '100:00:00';% max: 100h
            OPTS.nlab     = nlab;       % total: 4
            
        otherwise
            error('Avaliable queues are "clexpress", "clmedium", "cllong", "clfocean');
    end

    OPTS.nodes        = OPTS.nlab;
    OPTS.procspernode = 16;
    OPTS.mem_gb       = min(128,2 * OPTS.nlab * OPTS.procspernode);
    OPTS.path2source  = pwd;
    OPTS.path2mutils  = '/sfs/fs6/home-geomar/smomw219/MATLAB/Matlab-Libs/mutils-0.4-2_openmp';
    OPTS.path2tools   = '/sfs/fs6/home-geomar/smomw219/MATLAB/Matlab-Libs/tools';
    OPTS.matlab_version = version('-release');
    
else
    error(' Only config==nec and config==local have been defined yet.');
end

end