function write_restart_file(SETTINGS,VAR,COMM,istep,iplot,time,dt,irestart) %#ok<INUSL>
% Usage: write_restart_file(SETTINGS,VAR,COMM,istep,iplot,time,dt,irestart)
%
% Purpose: Load restart data of previous model calculation.
%
% Input:
%   SETTINGS : [structure] : model parameters
%   VAR      : [structure] : major variable fields
%   istep    : [scalar]    : current number of time step
%   iplot    : [scalar]    : current plot file number
%   time     : [scalar]    : current simulation time
%   dt       : [scalar]    : length of current time step
%   irestart : [scalar]    : number of restart file (1 or 2)
%
% Output:
%   none
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

fprintf(1, ' Writing restart file number #%1i...',irestart);
filename = [SETTINGS.outdir '/' COMM.prefix '_RestartFile' num2str(irestart)];
save(filename,'VAR','dt','time','istep','iplot','dt');
fprintf(' done.\n');

end % END OF FUNCTION write_restart_file