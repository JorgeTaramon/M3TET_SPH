function detJ = elblk_detA(Jx,Jy,Jz)

detJ = Jx(:,1).*Jy(:,2).*Jz(:,3) ...
     + Jx(:,2).*Jy(:,3).*Jz(:,1) ...
     + Jx(:,3).*Jy(:,1).*Jz(:,2) ...
     - Jx(:,3).*Jy(:,2).*Jz(:,1) ...
     - Jx(:,1).*Jy(:,3).*Jz(:,2) ...
     - Jx(:,2).*Jy(:,1).*Jz(:,3);
      
end