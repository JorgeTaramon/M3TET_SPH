function [Bx,By,Bz] = elblk_transpose_A(Ax,Ay,Az)

% Calculates B=A', where A is stored in MILAMIN-style block format
%
% Storage of 3x3 element matrix A in element-block format:
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
% Now we want to calculate B = [Bx(iel,:)' By(iel,:)' Bz(iel,:)']
% such that B = A';
% 
% Hence:
% Bx(iel,:) = [Ax(iel,1) Ay(iel,1) Az(iel,1)];
% By(iel,:) = [Ax(iel,2) Ay(iel,2) Az(iel,2)];
% Bz(iel,:) = [Ax(iel,3) Ay(iel,3) Az(iel,3)];

Bx = [Ax(:,1) Ay(:,1) Az(:,1)];
By = [Ax(:,2) Ay(:,2) Az(:,2)];
Bz = [Ax(:,3) Ay(:,3) Az(:,3)];

end