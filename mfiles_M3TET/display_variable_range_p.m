function display_variable_range_p(VAR,COMM,fidl)
% Usage: display_variable_range_p(VAR,COMM,fidl)
%
% Purpose: Shows min/max/mean values of all variables stored in structure
%          "VAR". Use cell-array "exclude" to skip variables.
%
% Input:
%   VAR    : [structure] : major variable fields, each is a vector
%   COMM   : [structure] : inter-subdomain communication data
%   fidl   : [scalar]    : file handle (optional)
%
% Output:
%         none
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH March 2013
%

if nargin<2
    COMM.nsd = 1;
end
if nargin<3
    fidl = 1;
end

exclude = {};

varnames = fieldnames(VAR);
nvar     = length(varnames);
fprintf(fidl,'\n Checking ranges of major variables...\n');
fprintf(fidl,' --------------------------------------------------------\n');
fprintf(fidl,' Variable  |      min     |       max     |      mean   |\n');
for i=1:nvar
    varname = varnames{i};

    if isstruct(VAR.(varname))
        continue
    end
    
    skip = 0;
    for j=1:length(exclude)
        if strcmp(varname,exclude{j})
            skip = 1; break
        end
    end
    
    if skip; continue; end
    
    nc = size(VAR.(varname),2);
    m  = 9 - length(varname);
    if nc>1
        m = m -1;
        for ic=1:nc
            fprintf(' %s%1i%s | ',varname,ic,repmat(' ',1,m));
            echo_minmaxmean(VAR.(varname)(:,ic),COMM,fidl);
        end
    else
        fprintf(fidl,' %s%s | ',varname,repmat(' ',1,m));
        echo_minmaxmean(VAR.(varname),COMM,fidl);
    end
end
fprintf(fidl,' --------------------------------------------------------\n\n');

end

% #########################################################################
%                               SUBFUNCTIONS
% #########################################################################

function echo_minmaxmean(v,COMM,fidl)

if COMM.nsd==1
    vmin  = min(v(:));
    vmax  = max(v(:));
    vmean = mean(v(:));
else
    vmin         = COMM.minLabs(v(:));
    vmax         = COMM.maxLabs(v(:));
    unique_nodes = COMM.unique_nodes{1}(COMM.unique_nodes{1}<=length(v));
    vmean        = COMM.sum_all(sum(v(unique_nodes))) ./ ...
                   COMM.sum_all(length(unique_nodes));
end
fprintf(fidl,'%+0.4e  |  %+0.4e  | %+0.4e |\n',vmin,vmax,vmean);
   
end
