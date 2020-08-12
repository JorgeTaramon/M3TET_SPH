function create_output_folder(SETTINGS)
% Usage: create_output_folder(SETTINGS)
%
% Purpose: Checks if output folder "SETTINGS.outdir" exists; 
%          creates it if not existing
%
% Input:
%   SETTINGS : [structure] : model parameters
%
% Output:
%   none
%
% JH Dec 2012
%

if labindex==1
    % IMPORTANT: Master worker must create output directory first,
    %            all other workers must wait
    try
        a = ls(SETTINGS.outdir);
    catch
        a = [];
    end
    if isempty(a)
        try
            mkdir(SETTINGS.outdir);
        catch
            error(' Could not create output directory. Check if parent directory exists.');
        end
    end
end

% Synchronize workers
labBarrier

end