function [detJ,invJx,invJy,invJz] = elblk_invJacobian_cartesian(GCOORD,EL2NOD,dN)

[nnodel,nel] = size(EL2NOD);

% COORDINATES OF CORNER/VERTEX NODES OF LL ELEMENTS IN BLOCK
ECOORD_x = reshape( GCOORD(1,EL2NOD), nnodel, nel );
ECOORD_y = reshape( GCOORD(2,EL2NOD), nnodel, nel );
ECOORD_z = reshape( GCOORD(3,EL2NOD), nnodel, nel );
Jx       = ECOORD_x'*dN;
Jy       = ECOORD_y'*dN;
Jz       = ECOORD_z'*dN;

detJ     = elblk_detA(Jx,Jy,Jz);

if any(detJ<0)
    iel = find(detJ<0,1,'first');

    figure(76);clf
    hold on
    for i=1:nnodel
        scatter3(ECOORD_x(i,iel),ECOORD_y(i,iel),ECOORD_z(i,iel),200,'k');
        text(ECOORD_x(i,iel),ECOORD_y(i,iel),ECOORD_z(i,iel),num2str(i),...
            'Fontsize',16);
    end
    view(3)
    if nnodel==4
        segs = [1 2; 2 3; 3 4; 1 4; 1 3];
    else
        segs = [1 5; 5 2;
                2 6; 6 3;
                3 7; 7 4;
                4 8; 8 1;
                1 9; 9 3;
                4 10; 10 2];
    end
    for i=1:size(segs,1)
        line(ECOORD_x(segs(i,:),iel),...
             ECOORD_y(segs(i,:),iel),...
             ECOORD_z(segs(i,:),iel),'Color','k');
    end
    view(3)

    error('negative Jacobi')
end

if nargout==4
    [invJx,invJy,invJz] = elblk_invert_A(Jx,Jy,Jz,detJ);
else
    invJx = [];
    invJy = [];
    invJz = [];
end

end