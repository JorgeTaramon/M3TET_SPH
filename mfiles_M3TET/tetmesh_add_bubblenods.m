function [GCOORD,EL2NOD] = tetmesh_add_bubblenods(GCOORD,EL2NOD,nbubble)

nnodel = size(EL2NOD,1);
nel    = size(EL2NOD,2);
nnod0  = max(EL2NOD(:));

switch nbubble
    case 1
        x_midnods = sum(reshape(GCOORD(1,EL2NOD(1:4,:)),4,nel))./4;
        y_midnods = sum(reshape(GCOORD(2,EL2NOD(1:4,:)),4,nel))./4;
        z_midnods = sum(reshape(GCOORD(3,EL2NOD(1:4,:)),4,nel))./4;
        GCOORD    = [GCOORD [x_midnods; y_midnods; z_midnods]];
        nnod      = size(GCOORD,2);
        EL2NOD(nnodel+1,:) = uint32(nnod0+1 : nnod);

    case 4
        error('To be coded.');
        
    otherwise
        error('Can only add 1 or 4 bubble nodes.');
end

end