function COMM = init_COMM()
% Usage: COMM = init_COMM()
% 
% Purpose: Initilizes structure COMM, which is needed in all parallel codes
%          (even when running serial).
%
% Input:
%   none
%
% Output:
%   COMM : [structure] : info on parallel (or serial) setup
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Jan 2012
% JH Feb 2016 : added new machines
%

COMM.nsd    = numlabs;  % number of subdomains running the code
COMM.myid   = labindex; % index of subdomain
COMM.prefix = ['Lab' num2str_d(COMM.nsd,2) 'x' num2str_d(COMM.myid,2)];

hostname = get_hostname;
if ismac
    switch hostname
        case 'Sneafell'
            COMM.ncores = 1;
        case 'jorgeMAC'
            COMM.ncores = 1; % Although MUTILS was compiled with use_openmp = 1
                             % ('mutils_config.m') MEX files are compiled
                             % with no Openmp (see README in mutils-0.4-2 folder)
%             [~,result]  = system('sysctl -n hw.ncpu');
%             COMM.ncores = str2double(result);
        otherwise
            % make sure MUTILS was compiled with with useopenmp = 1 (see
            % file 'mutils_config.m') and works!
            COMM.ncores = 1;
    end
    
elseif isunix
    % Names of unix machines on wwhich hyperthreading is enabled (operating
    % system 'simulating' twice the number of available physical CPUs; this
    % has performance benefits for smaller desktop activities but not for
    % computational expensive jobs)
    unixservers_with_hyperthread = {'nesh-fe'};
    % Check if thsi computer is in the above list
    hyperthread = 0;
    if ismember(hostname,unixservers_with_hyperthread)
        hyperthread = 1;
    end
    % Get number of cores on this computer
    [~,result]  = system('numactl --show | grep physcpubind');
    % result = 'physcpubind: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15'
    COMM.ncores = length(str2num(result(strfind(result,'0') : end))); %#ok<ST2NM>
    if hyperthread
        COMM.ncores = floor(COMM.ncores/2);
    end
    
elseif ispc
    switch hostname
        case 'b3pc9'
            COMM.ncores = 1;
        case 'b3pc14'
            COMM.ncores = 4;
        case 'b3pc16'
            COMM.ncores = 16;
        case 'b3pc17'
            COMM.ncores = 8;
        case 'Colima'
            COMM.ncores = 10;
        case 'Maupiti2'
            COMM.ncores = 1; % because MUTILS' sparse_convert is terribly slow on windows
        otherwise
            % make sure MUTILS was compiled with with useopenmp = 1 (see
            % file 'mutils_config.m') and works!
            COMM.ncores = 1;
    end
end

if numlabs==1
    COMM.master   = 1;
    COMM.slaves   = [];
    COMM.nthreads = COMM.ncores; % for MUTILS' spmv
    return % if nsd>1 but numlabs==1 we're in debugging mode. RETURN
end

SD_hosts = cell(COMM.nsd,1);
for ilab=1:COMM.nsd
    if COMM.myid==ilab
        SD_hosts{ilab} = labBroadcast( COMM.myid, hostname);
    else
        SD_hosts{ilab} = labBroadcast( ilab );
    end
end

labs_on_shared_mem = [];
for ilab=1:COMM.nsd
    if hostname==SD_hosts{ilab}
        labs_on_shared_mem = [labs_on_shared_mem ilab]; %#ok<AGROW>
    end
end

COMM.master   = min( labs_on_shared_mem );
COMM.slaves   = setdiff( labs_on_shared_mem , COMM.master );

COMM.nthreads = floor(COMM.ncores / length(labs_on_shared_mem)); % for MUTILS' spmv

end