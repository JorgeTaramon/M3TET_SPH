function SETTINGS = open_logfile(SETTINGS,COMM,on_screen)

% IMPORTANT: Master node must create output directory first,
% all other nodes must wait
labBarrier

% Name of log file
log_file = [SETTINGS.outdir filesep COMM.prefix '.log'];

if (numlabs>1 || ~isempty(getCurrentWorker)) && ~on_screen
    if SETTINGS.restart
        fidl = fopen(log_file,'a');
    else
        fidl = fopen(log_file,'w');
    end
    if fidl<0
        error(' Cannot create log file\n "%s"',log_file);
    end
else
    fidl = 1;
end
SETTINGS.fid_log = fidl;

fprintf(fidl,'\n\n\n\n');
fprintf(fidl,' =========================================================================\n');
if SETTINGS.restart
fprintf(fidl,'                       CONTINUING A CALCULATION \n');
else
fprintf(fidl,'                       STARTING NEW CALCULATION \n');
end
if COMM.nsd==1
    fprintf(fidl,' RUNNING IN SERIAL MODE. THIS IS WORKER "%s"\n\n',...
        get_hostname);
else
    fprintf(fidl,' RUNNING IN PARALLEL MODE. THIS IS WORKER %1i of %1i ("%s").\n',...
        COMM.myid,COMM.nsd,get_hostname);
end
fprintf(fidl,' DATA WILL BE WRITTEN TO FOLDER:\n    "%s"\n',SETTINGS.outdir);
% Write structure "COMM" to log file
fprintf(fidl,'\n SETTINGS SHARED MEMORY OPERATIONS (OpenMP)\n');
fprintf(fidl,' ncores        : %3i\n',COMM.ncores);
fprintf(fidl,' nthreads      : %3i\n',COMM.nthreads);
fprintf(fidl,' master worker : %3i\n',COMM.master);
fprintf(fidl,' slave workers :');
for i=1:length(COMM.slaves)
    fprintf(fidl,' %3i',COMM.slaves(i));
    if mod(i,14)==0
        fprintf(fidl,'\n                ');
    end
end
fprintf(fidl,'\n');
fprintf(fidl,' =========================================================================\n');

end