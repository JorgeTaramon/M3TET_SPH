function jm = configJM(jm,OPTS)

% Check input arguments
% #####################
check_options(nargin,OPTS); % *SUBFUNCTION*

% Construct input string for jobmanager
% #####################################
if strcmp(OPTS.config,'rzcl')
    % jm.ResourceTemplate (examples)
    % local    : '-l select=12';
    % rzcluster: '-l select=32:mem=2gb -l walltime=48:00:00 -q small -oe';
    %
    % jm.SubmitArguments (examples)
    % local    : '-oe';
    % rzcluster: '-q small -oe';

    ResourceTemplate = sprintf('-l select=%1i',OPTS.nlab);
    if isfield(OPTS,'mem_gb')
        ResourceTemplate = sprintf('%s:mem=%1igb',ResourceTemplate,OPTS.mem_gb);
    end
    if isfield(OPTS,'walltime')
        ResourceTemplate = sprintf('%s -l walltime=%s',ResourceTemplate,OPTS.walltime);
    end
    jm.ResourceTemplate = ResourceTemplate;

    if isfield(OPTS,'queue')
        SubmitArguments = sprintf('-q %s ',OPTS.queue);
    else
        SubmitArguments = '';
    end
    if isfield(OPTS,'SubmitArguments')
        SubmitArguments = sprintf('%s%s',SubmitArguments,OPTS.SubmitArguments);
    end
    jm.SubmitArguments = SubmitArguments;

elseif strcmp(OPTS.config,'nec')
    % Some additional fields are now generated based on the available
    % fields in "OPTS". Everything needed will be stored in structure
    % "PBS", which is then written to file "start_pbs.m". This is currently
    % required to start jobs on the NEC cluster.
    
    % Save OPTS variables to script
    % Note that a) all PBS.* variables need to be defined for proper execution, and
    %           b) the filename 'start_pbs.m' is hard-coded sofar.
    PBS.path         = OPTS.path2source;
    PBS.config       = 'generic';
    PBS.queue        = OPTS.queue;
    PBS.nodes        = OPTS.nodes;
    PBS.procspernode = OPTS.procspernode;
    PBS.elapstime    = OPTS.walltime;

    ind_col = find(OPTS.walltime==':');
    h       = str2double(OPTS.walltime(1:ind_col(1)-1));
    m       = str2double(OPTS.walltime(ind_col(1)+1:ind_col(2)-1));
    s       = str2double(OPTS.walltime(ind_col(2)+1:end));
    s_cpu   = OPTS.nodes * (3600*h + 60*m + s);
    h_cpu   = floor(s_cpu/3600);
    s_cpu   = s_cpu-h_cpu*3600;
    m_cpu   = floor(s_cpu/60);
    s_cpu   = s_cpu-m_cpu*60;
    cputime = sprintf('%2i:%2i:%2i',h_cpu,m_cpu,s_cpu);
    ind0    = cputime==' ';
    cputime(ind0) = '0';
    PBS.cputime = cputime; % accumulated CPU time per node

    PBS.mem     = [num2str(OPTS.mem_gb) 'gb']; %#ok<STRNU>
    matlab.io.saveVariablesToScript('start_pbs.m','PBS');
    
end

end % END OF FUNCTION configJM

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function check_options(n,OPTS)

if n<2
    error('Function must be called "configJM(jm,OPTS)"');
end

if ~isstruct(OPTS)
    error('Configuration must be defined in a structure: OPTS.config, OPTS.nlab,...');
end

if ~isfield(OPTS,'config')
    error('OPTS.config must be the name of the configuration.');
end

if ~isfield(OPTS,'nlab')
    error('OPTS.nlab must be the number of Matlab workers for the job');
end

switch OPTS.config
case 'rzcl'
    error('config==rzcl need to be checked again.');
    if ~isfield(OPTS,'mem_gb')
        error('OPTS.mem_gb must be defined. E.g.: "OPTS.mem_gb = 2"');
    end
    if ~isfield(OPTS,'walltime')
        error('OPTS.walltime must be defined. E.g.: "OPTS.walltime = 48:30:00"');
    elseif isnumeric(OPTS.walltime)
        error('OPTS.walltime must be defined as a string. E.g.: "OPTS.walltime = 48:30:00"');
    end
    if ~isfield(OPTS,'queue')
        error('OPTS.queue must be defined. Avaliable queues are "foexpress", "f_ocean2"');
    else
        if strcmp(OPTS.queue,'f_ocean2')
            disp(' Using queue "f_ocean2"');
        elseif strcmp(OPTS.queue,'foexpress')
            disp(' Using queue "f_ocean2"');
        else
            error('OPTS.queue must be "foexpress" or "f_ocean2"');
        end
    end

case 'nec'
    if ~isfield(OPTS,'mem_gb') || ~isnumeric(OPTS.mem_gb)
        error('"OPTS.mem_gb" must be defined as a number. E.g.: "OPTS.mem_gb = 64"');
    end
    
    if ~isfield(OPTS,'nodes') || ~isnumeric(OPTS.nodes)
        error('"OPTS.nodes" must specify the number of nodes (each has 16 cores!)');
    end
    
    if ~isfield(OPTS,'procspernode') || ~isnumeric(OPTS.procspernode)
        txt = ['"OPTS.procspernode" must specify the number of processes per node.\n' ...
               'Note that each node has 16 cores and total number requested cores will be \n' ...
               'OPTS.nodes * OPTS.procspernode!'];
        error(txt);
    end
    
    if OPTS.nodes * OPTS.procspernode > 128
        error('OPTS.nodes * OPTS.procspernode > 128 (max number of licenses)!');
    end
    
    if ~isfield(OPTS,'queue') || ~ismember(OPTS.queue,{'clexpress','clmedium','cllong'})
        error('OPTS.queue must be defined. Avaliable queues are "clexpress", "clmedium", "cllong"');
    else
        fprintf('\n Using queue "%s". Current usage:\n',OPTS.queue);
        cmd        = sprintf('qstatall -q %s',OPTS.queue);
        [~,result] = system(cmd);
        disp(result); fprintf('\n');
    end

    if ~isfield(OPTS,'walltime')
        error('OPTS.walltime must be defined. Format: "OPTS.walltime = hh:mm:ss"');
    elseif isnumeric(OPTS.walltime)
        error('OPTS.walltime must be defined as a string. Format: "OPTS.walltime = hh:mm:ss"');
    else
        ind_col = find(OPTS.walltime==':');
        if length(ind_col)~=2
            error('OPTS.walltime seems to have the wrong format. Format: "OPTS.walltime = hh:mm:ss"')
        end
    end
end

end % END OF SUBFUNCTION check_options