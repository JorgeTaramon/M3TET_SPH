function START_PARALLEL_JOB(OPTS,mainfile,OPTS_RUN)
% Usage: START_PARALLEL_JOB(OPTS,mainfile,OPTS_RUN)
% 
% Purpose: Find and configure scheduler/jobmanager and submit parallel job.
%
% Input:
%   OPTS     : [structure] : must contain number of workers & configuration
%   mainfile : [string]    : name of file to be executed by workers
%   OPTS_RUN : [structure] : input argument for main file
%
% Output:
%   none
%
% Content of structure "OPTS"
% ...fields that are always required:
%      OPTS.config      % jobmanager configuration (see examples below)
%      OPTS.nlab        % number of workers (==number of subdomains)
%      OPTS.path2source % path to function "mainfile"
%      OPTS.path2mutils % path to mutils installation
%      OPTS.path2tools  % path to folder with supporting m-files
%
% ...optional fields
%      OPTS.nthread  = <number of threads>
%        % number of threads that each worker is allowed to use for OpenMP
%        % example: nlab=4, nthread=2 --> 4 subdomains, each 2 threads
%                   --> machine should have >=8 free cores
%
% Note that there is a function "DEFINE_OPTS" to quickly create "OPTS" with
% some standard settings that can be edited.
%
% Example for LOCAL scheduler:
%      OPTS.config   = 'local';
%      OPTS.nlab     = 4;
%
% Example for RZCLUSTER:
%      OPTS.config   = 'rzcl';
%      OPTS.nlab     = 16;         % 1-128
%      OPTS.mem_gb   = 2;          % 
%      OPTS.walltime = '80:00:00'; % 
%      OPTS.queue    = 'medium';   % 'express' 'small' 'medium' 'long'
%
% Example for NEC:
%      OPTS.config   = 'nec';
%      OPTS.nlab     = 8;           % 1-128
%      OPTS.nnodes   = 2;           % number of nodes (each has 16 cores)
%      OPTS.cores_per_node = 16;    % number of cores
%         note that OPTS.nnodes*OPTS.cores_per_node != OPTS.nlab
%      OPTS.mem      = '64gb';      % required memory
%      OPTS.walltime = '80:00:00';  % express:2h, medium:48h, long:100h
%      OPTS.queue    = 'clexpress'; % 'clexpress', 'clmedium', 'cllong'
%      
%
% Joerg Hasenclever, jhasenclever@geomar.de
%
% JH Jul 2013
% JH Nov 2014 : cleaned up
% JH Dec 2014 : removed abakus; added NEC-cluster
% JH Jan 2016 : added option to set number of threads; improved matlab
%               version check
%

% If this function is called by the workers (see subfunction), the 
% actual "mainfile" with optional input argument will be executed.
% ###################################################################
if nargin==1 && isfield(OPTS,'cmd4workers')
    run_parallel_cmd(OPTS); % *SUBFUNCTION*
    return
end

% Check input
% ###########
if nargin<1
    error('Function must be called "START_PARALLEL_JOB(OPTS,mainfile)"');
end
if nargin==3 && ~isstruct(OPTS_RUN)
    error('Function must be called "START_PARALLEL_JOB(OPTS,mainfile,OPTS_RUN)", where OPTS_RUN is a structure.');
end
OPTS = check_input(OPTS,nargin); % *SUBFUNCTION*
if isempty(ls([OPTS.path2source '/' mainfile '.m']))
    error('Main file %s.m not found.',mainfile);
end
if nargin<3 && OPTS.no_questions==0
    s = input(sprintf(' No input argument for %s has been defined. Continue (Y/N) ?',mainfile),'s');
    if s~='y' && s~='Y'; error('Aborted by user.'); end
end

% Find out which Matlab version we're using
% #########################################
OPTS.matlab_version    = version('-release');
OPTS.matlab_release_yr = str2double(OPTS.matlab_version(1:4));

% Find jobmanager
% ###############
jm = findJM(OPTS);

% Configure jobmanager
% ####################
switch OPTS.config
case 'local'
    s = input(sprintf(' Will start job with %1i workers using local jobmanager.\n Continue (Y/N) ?',OPTS.nlab),'s');
case 'rzcl'
    jm = configJM(jm,OPTS);
    fprintf(' Will submit the following options to JM on rzcluster: %s %s\n',...
        jm.ResourceTemplate,jm.SubmitArguments);
    s = input(' Please verify the above submit arguments.\n Continue (Y/N)? \n Answer: ','s');
case 'nec'
    jm = configJM(jm,OPTS);
    fprintf('\n The following configuration was written to file "start_pbs.m":\n\n')
    system('more start_pbs.m');
    s = input(' Please verify the above submit arguments.\n Continue (Y/N)? \n Answer: ','s');
end
if ~(strcmp(s,'y') || strcmp(s,'Y'))
    error(' Aborted by user.');
end

% Command to be executed by each Matlab worker
% ############################################
SUBMISSION.path2source = OPTS.path2source;
if nargin==2
    SUBMISSION.cmd4workers = sprintf('%s;',mainfile);
elseif nargin==3
    SUBMISSION.OPTS_RUN    = OPTS_RUN;
    SUBMISSION.cmd4workers = sprintf('%s(OPTS_RUN);',mainfile);
end
fprintf(' Each worker will execute the following command:\n       "%s"\n',...
    SUBMISSION.cmd4workers);
if OPTS.no_questions==0
    s = input(' Continue (Y/N)? ','s');
    if ~(strcmp(s,'y') || strcmp(s,'Y'))
        error(' Aborted by user.');
    end
end

% Create parallel job
% ###################
paths = {SUBMISSION.path2source};
if ~isempty(OPTS.path2mutils)
    paths        = addpaths_mutils(paths,OPTS.path2mutils); % *SUBFUNCTION*
end
if ~isempty(OPTS.path2tools)
    paths{end+1} = OPTS.path2tools;
end

if OPTS.matlab_release_yr>2013
    job                 = createCommunicatingJob(jm,'Type','spmd');
    job.NumWorkersRange = [OPTS.nlab OPTS.nlab];
    job.AdditionalPaths = paths;
else
    job = createParallelJob(jm,'MinimumNumberOfWorkers',OPTS.nlab,...
                               'MaximumNumberOfWorkers',OPTS.nlab,...
                               'PathDependencies',paths); %#ok<DCRTPJ>
end


% Create the task(s) of the parallel job
% ######################################
task = createTask(job, @START_PARALLEL_JOB, 0, {SUBMISSION}); %#ok<NASGU>
% syntax:
% task = createTask(job, F, N, {inputargs})
% task : Task object or vector of task objects.
% job  : The job that the task object is created in.
% F    : A handle to the function that is called when the task is
%        evaluated, or an array of function handles.
% N    : The number of output arguments to be returned from execution of
%        the task function. This is a double or array of doubles.
% {inputargs} : A row cell array specifying the input arguments to be 
%               passed to the function F. Each element in the cell array
%               will be passed as a separate input argument.


% Submit the job to the workers
% #############################
fprintf(' Submitting job...');
submit(job);
fprintf('done.\n');
if strcmp(OPTS.config,'nec')
    cmd        = sprintf('qstatall -q %s',OPTS.queue);
    [~,result] = system(cmd);
    disp(result);
end
fprintf('\n -------------------------------\n');

end % END OF FUNCTION START_PARALLEL_JOB

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function OPTS = check_input(OPTS,n)

if ~isfield(OPTS,'config')
    error('OPTS.config (type of cluster profile) is not defined.');
end
if ~isfield(OPTS,'nlab')
    error('OPTS.nlab (number of workers) is not defined.');
end
if ~isfield(OPTS,'path2source')
    error('OPTS.path2source (path to function to be executed) is not defined.');
end
if ~isfield(OPTS,'path2mutils')
    error('OPTS.path2mutils must be defined (or set empty).');
end
if ~isfield(OPTS,'path2tools')
    error('OPTS.path2tools must be defined (or set empty).');
end
if n<2
    error('Name of main file to be executed by workers must be 2nd input argument');
end
if ~isfield(OPTS,'no_questions')
    OPTS.no_questions = 0;
end

end % END OF SUBFUNCTION check_input

% #########################################################################

function run_parallel_cmd(SUBMISSION)

if ~isfield(SUBMISSION,'path2source') || isempty(SUBMISSION.path2source) || ...
   ~ischar(SUBMISSION.path2source)
    error('"SUBMISSION.path2source" is not defined correctly.');
end
if ~isfield(SUBMISSION,'cmd4workers') || isempty(SUBMISSION.cmd4workers) || ...
   ~ischar(SUBMISSION.cmd4workers)
    error('"SUBMISSION.cmd4workers" is not defined correctly.');
end

% Make all workers "cd" into the folder where the "mainfile" is located
cd(SUBMISSION.path2source);
labBarrier

% Let all workers start the calculation by calling "mainfile"
if isfield(SUBMISSION,'OPTS_RUN')
    OPTS_RUN = SUBMISSION.OPTS_RUN; %#ok<NASGU>
end
eval(SUBMISSION.cmd4workers);

end % END OF SUBFUNCTION run_parallel_cmd

% #########################################################################

function paths = addpaths_mutils(paths,path2mutils)

paths{end+1} = path2mutils;
paths{end+1} = [path2mutils '/triangle'];
paths{end+1} = [path2mutils '/SuiteSparse'];
paths{end+1} = [path2mutils '/mutils'];
paths{end+1} = [path2mutils '/mutils/quadtree'];
paths{end+1} = [path2mutils '/mutils/libutils'];
paths{end+1} = [path2mutils '/mutils/libmatlab'];
paths{end+1} = [path2mutils '/mutils/interp'];
paths{end+1} = [path2mutils '/mutils/reorder'];
paths{end+1} = [path2mutils '/mutils/sparse'];

end