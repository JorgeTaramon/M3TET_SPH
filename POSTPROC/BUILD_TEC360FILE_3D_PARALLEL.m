function BUILD_TEC360FILE_3D_PARALLEL(nlab,OPTS)
% Usage: BUILD_TEC360FILE_3D_PARALLEL(nlab,OPTS)
% 
% Purpose: Find and configure local scheduler and submit parallel job.
%
% Input:
%   nlab     : [integer]   : number of CPUs to process job
%   OPTS     : [structure] : must contain number of workers & configuration
%
% Output:
%   none
%
% Part of HT3_SP - 3D FINITE ELEMENT SINGLE PHASE POROUS FLOW CODE
%
% JH Jul 2013
%


% Check input
% ###########
if nargin<1
    error('Function must be called "START_PARALLEL_JOB(nlab)"');
end
if nargin<2
    OPTS = [];
end
if ~isfield(OPTS,'wait_until_finished')
    OPTS.wait_until_finished = 0;
end

% Find local scheduler
% ####################
jm = findResource('scheduler','type','local');
jm.ClusterSize =  12;

% % Configure scheduler/jobmanager
% % ##############################
% s = input(sprintf(' Will start job on %1i workers using local scheduler. Continue (Y/N) ?',nlab),'s');
% if ~(strcmp(s,'y') || strcmp(s,'Y'))
%     error(' Aborted by user.');
% end

% Create parallel job
% ###################
job   = createParallelJob(jm,'MinimumNumberOfWorkers',nlab,...
                             'MaximumNumberOfWorkers',nlab,...
                             'FileDependencies',{'BUILD_TEC360FILE_3D'});


% Create the task(s) of the parallel job
% ######################################
% Input is given to worker by function arguments    
task = createTask(job, @BUILD_TEC360FILE_3D, 0, {OPTS}); %#ok<NASGU>
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
submit(job);
disp(' Job submitted, use these command to view results')
disp(' jm     = findJM("local")');
disp(' job    = findJob(jm)')
disp(' task   = findTask(job)');
disp(' output = getAllOutputArguments(job)');

if OPTS.wait_until_finished
    while 1
        pause(30);
        state = get(job, 'State');
        if strcmp(state,'finished');
            break
        end
    end
end

end