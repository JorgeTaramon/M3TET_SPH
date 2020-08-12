function [A_uniq,ia2,ja2] = unique_keep_order(A)
% Usage: [A_uniq,ia2,ja2] = unique_keep_order(A)
%
% Purpose: Creates a list of unique rows in matrix A without changing the
%          order of the rows.
%
% Input:
%   A       : [matrix] : matrix of which the unique rows shall be returned
%
% Output:
%   A_uniq  : [matrix] : matrix with the unique rows of A, without
%                        re-ordering the rows
%   ia2     : [vector] : pointer A_uniq --> A
%   ja2     : [vector] : pointer A --> A_uniq
%
% Written by J.Hasenclever, 2007-2010
% Email contact: jhasenclever@geomar.de
%
% JH Mar 2016
%

% find FIRST occurrence of unique rows in A
% use sort so that values in each row are sorted increasingly
[~,ia,ja]   = unique(sort(A,2),'rows','first');

% create index ia2 to extract edges_uniq from A
[ia2,~,ja2] = unique(ia(ja));

% extract unique edges from input
A_uniq  = A(ia2,:);

end % END OF FUNCTION unique_keep_order