function write_output_file(SETTINGS,VAR,COMM,time,istep,dt,iplot)  %#ok<INUSL>
% Usage: write_output_file(SETTINGS,VAR,COMM,time,istep,dt,iplot)
% 
% Purpose: Saves major variables as mat-file in output folder.
%
% Input:
%   SETTINGS : [structure] : model parameters
%   VAR      : [structure] : major variable fields
%   time     : [scalar]    : current simulation time
%   istep    : [scalar]    : current number of time step
%   dt       : [scalar]    : length of current time step
%   iplot    : [scalar]    : current plot file number
%
% Output:
%   none (saves Matlab data to output folder)
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Dec 2012
%

fid_log = SETTINGS.fid_log;
fprintf(fid_log,' Writing output file #%1i...',iplot);
filename = [SETTINGS.outdir '/' COMM.prefix '_Data_' num2str_d(iplot,4) '.mat'];
save(filename,'istep','time','dt','VAR');
fprintf(fid_log,'done.\n');

end % END OF FUNCTION write_output_file