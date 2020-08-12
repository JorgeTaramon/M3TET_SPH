function [Cx,Cy,Cz] = elblk_AxB(Ax,Ay,Az,Bx,By,Bz)

% Calculates C = A * B, where A, B, and C are in MILAMIN-style block format

% Storage of 3x3 element matrices A and B in element-block format:
%   Ax(1,:) is A(:,1) for the first element in the block (here called "iel")
%   Ay(1,:) is A(:,2) for the first element in the block
%   Az(1,:) is A(:,3) for the first element in the block
%
% Hence for iel: A = [Ax(iel,:)' Ay(iel,:)' Az(iel,:)']
%                     _                             _
%                    | Ax(iel,1) Ay(iel,1) Az(iel,1) |
%                  = | Ax(iel,2) Ay(iel,2) Az(iel,2) |
%                    |_Ax(iel,3) Ay(iel,3) Az(iel,3)_|
% 
% If matrix B is stored accodingly, then C = A * B is equal to
%    _                             _     _                             _
%   | Ax(iel,1) Ay(iel,1) Az(iel,1) |   | Bx(iel,1) By(iel,1) Bz(iel,1) |
%   | Ax(iel,2) Ay(iel,2) Az(iel,2) | * | Bx(iel,2) By(iel,2) Bz(iel,2) |
%   |_Ax(iel,3) Ay(iel,3) Az(iel,3)_|   |_Bx(iel,3) By(iel,3) Bz(iel,3)_|
% 
% If Cx(iel,:) is the frst column of the matrix C for element "iel" then:
% Cx(iel,:) = [Ax(iel,1)*Bx(iel,1)+Ay(iel,1)*Bx(iel,2)+Az(iel,1)*Bx(iel,3) ...
%              Ax(iel,2)*Bx(iel,1)+Ay(iel,2)*Bx(iel,2)+Az(iel,2)*Bx(iel,3) ...
%              Ax(iel,3)*Bx(iel,1)+Ay(iel,3)*Bx(iel,2)+Az(iel,3)*Bx(iel,3)]
% Cy(iel,:) is the same but using By instead of Bx
% and
% Cz(iel,:) is the same but Bz instead of Bx.

if nargin==0
    A   = [1 2 3; 4 5 6; 7 8 9];
    B   = 100 + A;
    C   = A*B;
    
    nn  = 4;
    Ax  = repmat(A(:,1)',nn,1);
    Ay  = repmat(A(:,2)',nn,1);
    Az  = repmat(A(:,3)',nn,1);
    Bx  = repmat(B(:,1)',nn,1);
    By  = repmat(B(:,2)',nn,1);
    Bz  = repmat(B(:,3)',nn,1);
end

Cx = [Ax(:,1).*Bx(:,1)+Ay(:,1).*Bx(:,2)+Az(:,1).*Bx(:,3) ...
      Ax(:,2).*Bx(:,1)+Ay(:,2).*Bx(:,2)+Az(:,2).*Bx(:,3) ...
      Ax(:,3).*Bx(:,1)+Ay(:,3).*Bx(:,2)+Az(:,3).*Bx(:,3)];
  
Cy = [Ax(:,1).*By(:,1)+Ay(:,1).*By(:,2)+Az(:,1).*By(:,3) ...
      Ax(:,2).*By(:,1)+Ay(:,2).*By(:,2)+Az(:,2).*By(:,3) ...
      Ax(:,3).*By(:,1)+Ay(:,3).*By(:,2)+Az(:,3).*By(:,3)];
  
Cz = [Ax(:,1).*Bz(:,1)+Ay(:,1).*Bz(:,2)+Az(:,1).*Bz(:,3) ...
      Ax(:,2).*Bz(:,1)+Ay(:,2).*Bz(:,2)+Az(:,2).*Bz(:,3) ...
      Ax(:,3).*Bz(:,1)+Ay(:,3).*Bz(:,2)+Az(:,3).*Bz(:,3)];

if nargin==0
    C - [Cx(1,:)' Cy(1,:)' Cz(1,:)']
    C - [Cx(nn,:)' Cy(nn,:)' Cz(nn,:)']
end

end