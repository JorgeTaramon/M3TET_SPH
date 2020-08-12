function [detJ,invJx,invJy,invJz] = calc_jacobian(GCOORD,EL2NOD,dN)

[nnodel,nel] = size(EL2NOD);

% COORDINATES OF CORNER/VERTEX NODES OF LL ELEMENTS IN BLOCK
ECOORD_x = reshape( GCOORD(1,EL2NOD), nnodel, nel );
ECOORD_y = reshape( GCOORD(2,EL2NOD), nnodel, nel );
ECOORD_z = reshape( GCOORD(3,EL2NOD), nnodel, nel );
Jx       = ECOORD_x'*dN;
Jy       = ECOORD_y'*dN;
Jz       = ECOORD_z'*dN;
detJ     = Jx(:,1).*Jy(:,2).*Jz(:,3) ...
         + Jx(:,2).*Jy(:,3).*Jz(:,1) ...
         + Jx(:,3).*Jy(:,1).*Jz(:,2) ...
         - Jx(:,3).*Jy(:,2).*Jz(:,1) ...
         - Jx(:,1).*Jy(:,3).*Jz(:,2) ...
         - Jx(:,2).*Jy(:,1).*Jz(:,3);

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
    % 1/det(Js)
    invdetJ    = 1./detJ;

    % Inversion of Jacobi matrices of all elements in block
    invJx      = zeros(nel,3);
    invJy      = zeros(nel,3);
    invJz      = zeros(nel,3);
    invJx(:,1) = invdetJ .* ( Jy(:,2).*Jz(:,3) - Jy(:,3).*Jz(:,2) );
    invJx(:,2) = invdetJ .* ( Jy(:,3).*Jz(:,1) - Jy(:,1).*Jz(:,3) );
    invJx(:,3) = invdetJ .* ( Jy(:,1).*Jz(:,2) - Jy(:,2).*Jz(:,1) );
    invJy(:,1) = invdetJ .* ( Jx(:,3).*Jz(:,2) - Jx(:,2).*Jz(:,3) );
    invJy(:,2) = invdetJ .* ( Jx(:,1).*Jz(:,3) - Jx(:,3).*Jz(:,1) );
    invJy(:,3) = invdetJ .* ( Jx(:,2).*Jz(:,1) - Jx(:,1).*Jz(:,2) );
    invJz(:,1) = invdetJ .* ( Jx(:,2).*Jy(:,3) - Jx(:,3).*Jy(:,2) );
    invJz(:,2) = invdetJ .* ( Jx(:,3).*Jy(:,1) - Jx(:,1).*Jy(:,3) );
    invJz(:,3) = invdetJ .* ( Jx(:,1).*Jy(:,2) - Jx(:,2).*Jy(:,1) );
else
    invJx = [];
    invJy = [];
    invJz = [];
end

end