function jm = findJM(OPTS)

% Find scheduler or jobmanager
% ############################
config = OPTS.config;
if OPTS.matlab_release_yr>2013
    if strcmp(config,'local')
        fprintf(' Connecting to Matlabs local jobmanager...');
        jm = parcluster('local');
        fprintf('done.\n');
        if OPTS.nlab>jm.NumWorkers
            fprintf(' This system has %1i (virtual) cores but you are requesting %1i workers.\n',...
                jm.NumWorkers,OPTS.nlab);
            fprintf(' Performance will be very poor!\n');
            fprintf(' ONLY CONTINUE if you are TESTING!\n');
            fprintf(' ONLY CONTINUE if NO OTHER JOBS run on this system!\n');
            s = input(' Continue (Y/N) ?','s');
            if ~(strcmp(s,'y') || strcmp(s,'Y'))
                error(' Aborted by user.');
            else
                jm.NumWorkers = OPTS.nlab;
            end
        end
        
    elseif strcmp(config,'nec')
        fprintf(' Connecting to NEC Jobmanager on NEC-cluster (RZ Uni Kiel)...');
        jm = parcluster('generic');
        fprintf('done.\n');
        jm.NumWorkers = OPTS.nlab;
        
    elseif strcmp(config,'rzcl')
        fprintf(' Connecting to PBSpro Jobmanager on rzcluster (RZ Uni Kiel)...');
        jm = parcluster('rzcl');
        fprintf('done.\n');
        jm.NumWorkers = OPTS.nlab;

    else
        error('Unknown configuration.')
    end
    
else % old commands for pre-2014 matlab
    if strcmp(config,'local')
        disp(' Using local scheduler');
        jm = findResource('scheduler','type','local');
        jm.ClusterSize =  12;

    elseif strcmp(config,'nec')
        error('Matlab R2014a must be used on NEC cluster.');

    elseif strcmp(config,'rzcl')
        disp(' Using PBSpro Jobmanager on rzcluster (RZ Uni Kiel)');
        jm = findResource('scheduler', 'type', 'pbspro',...
            'Configuration','rzcl');

    else
        error('Unknown configuration.')
    end
end

end % END OF FUNCTION findJM