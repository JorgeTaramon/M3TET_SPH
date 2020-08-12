function [invJx,invJy,invJz] = elblk_invert_A(Jx,Jy,Jz,detJ)

nel        = length(detJ);

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

end