function write_tecplot_data_3d(MESH,VAR,SETTINGS,NUMSCALE,time,istep,iplot,varnames)

tecfile  = [SETTINGS.outdir '/SphereData_' num2str_d(iplot,4) '.plt'];
nvar     = length(varnames);
EL2NOD   = MESH.EL2NOD{1};
nnodel   = size(EL2NOD,1);
if nnodel>4
    EL2NOD_split = tetmesh_p2_to_p1([],EL2NOD);
else
    EL2NOD_split = EL2NOD;
end
nnod        = max(EL2NOD_split(:));
GCOORD      = MESH.GCOORD(:,1:nnod);

% title of the plot, string
tdata.title = 'M3TET_SPH Data';

% integer, storing number of variables
tdata.Nvar  = nvar+3;

% cell array of variable names
tdata.varnames = [{'X' 'Y' 'Z'} varnames(:)'];

% integer vector array storing dataformat for each variable (1:Nvar).
% 1 for float; 2 for double; 3 for longInt; 4 for shortInt; 5 for Byte; 6 for Bit
tdata.vformat = 1*ones(1,3+nvar);

% tdata.FEvolumes % an array of finite element 3D volume data with 
%                   tetrahedron elements (4 points per element) or
%                   brick elements (8 points per element) defined on 3D
%                   XYZ coordinates. Each finite element volume is
%                   treated as a zone. FEvolumes contain the
%                   following information to define the zones:

% % strandID of the zone. Optional, default = -2
% tdata.FEvolumes.strandID = int32(iplot);

% zone name of each line
tdata.FEvolumes.zonename = sprintf('Plot %1i, Step %1i, Time %8.4f Myr',...
    iplot,istep,time);

% solution time associated with the zone Optional, default =0
tdata.FEvolumes.solutiontime = time;

% element connectivity list
tdata.FEvolumes.e2n = EL2NOD_split';

% x coordinate (1D array)
tdata.FEvolumes.x = GCOORD(1,:);

% y coordinate (1D array)
tdata.FEvolumes.y = GCOORD(2,:);

% z coordinate (1D array)
tdata.FEvolumes.z = GCOORD(3,:);

% is it block or point (deprecated and set to zero by default)
tdata.FEvolumes.datapacking = 1;

% location of variable: 0 --> nodal, 1 --> cell center (center of element))
tdata.FEvolumes.varloc = 0;

% value of each variable (nodal or cell)
% [nv,Nne], where nv is number of variables in v, Nne=Nn if nodal.
% Nne=Ne if cell-center Nn is number of nodes Ne is number of elements
for ivar=1:nvar
    varname = varnames{ivar};
    data    = VAR.(varname);
    nval    = length(data);
    if nval==MESH.nVnod
        data = interp_edge_node_values(EL2NOD,data);
    end
    switch varname
        case 'P'
            data = data .* NUMSCALE.P0 .* 1e-9;

        case 'Visc'
            data = max(data,1e-10); % cut-off to allow  using log10
            data = log10(data.*NUMSCALE.Visc0);
    end
    
    tdata.FEvolumes.v(ivar,:) = data(:)';
end


mat2tecplot(tdata,tecfile);

% WRITE A LIST OF ALL CREATED TECPLOT FILES (USEFUL TO COPY & PASTE INTO A
% LAYOUT FILE)
if iplot==0
    fid = fopen([SETTINGS.outdir '/tecplot_file_list.txt'],'w');
    s   = '"X" "Y" "Z"';
    for ivar=1:nvar
        s = [s ' "' varnames{ivar} '"'];
    end
    fprintf(fid,[s '\n']);
else
    fid = fopen([SETTINGS.outdir '/tecplot_file_list.txt'],'a');
end
outdir = SETTINGS.outdir;
outdir(outdir=='/') = '\';
if outdir(end)~='\'; outdir(end+1)='\'; end
ind     = strfind(outdir,'\');
i1      = ind(end-1)+1;
i2      = ind(end);
fwrite(fid,['"' outdir(i1:i2) 'SphereData_' num2str_d(iplot,4) '.plt" ']);
fclose(fid);

end % END OF FUNCTION write_tecplot_data_3d

% #########################################################################
%                               SUBFUNCTIONS
% #########################################################################

function data_all = interp_edge_node_values(EL2NOD,data_vert)

% data_vert : variable at vertex nodes
% data_edge : variable interpolated at edge nodes
% data_all  : variable at both edge and vertex nodes
           
nnod      = max(EL2NOD(:));
nVnod     = length(data_vert);
data_all  = zeros(nnod,1);
data_all(1:nVnod)        = data_vert;
data_all(EL2NOD(5:10,:)) = 0.5*(  data_vert(EL2NOD([1 2 3 4 1 2],:)) ...
                                + data_vert(EL2NOD([2 3 4 1 3 4],:)));

end