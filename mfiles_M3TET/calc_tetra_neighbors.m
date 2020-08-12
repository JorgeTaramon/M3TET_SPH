function NEIGHBORS = calc_tetra_neighbors(EL2NOD,iel)

nel          = size(EL2NOD,2);
if nargin==1
    iel = 1:nel;
else
    nel = length(iel);
end
% make pointer: face --> vertex nodes
FACE2NOD     = reshape(EL2NOD([2 3 4; 1 3 4;1 2 4;1 2 3]',iel),3,[])';
[~,~,eF2gF]  = unique(sort(FACE2NOD,2),'rows'); %#ok<UDIM>
clear FACE2NOD
% make pointer: element face --> global face
eF2gF        = reshape(eF2gF,4,nel);
% make pointer: global face --> elements
sub          = eF2gF(:);
val          = ((1:nel)' * ones(1,4))';
val          = val(:);
FACE2EL(1,sub) = val;
sub          = sub(end:-1:1);
val          = val(end:-1:1);
FACE2EL(2,sub) = val;
% now make pointer: element --> neighbors
NEIGHBORS    = zeros(4,nel,'uint32');
for iface=1:4
    tmp = FACE2EL(:,eF2gF(iface,:));
    ind = tmp~=[1:nel ; 1:nel];
    el  = find(any(ind>0,1));
    nb  = tmp(ind);
    ind = sub2ind([4 nel],repmat(iface,1,length(el)),el);
    NEIGHBORS(ind) = nb;
end

end % END OF SUBFUNCTION calc_tetra_neighbors