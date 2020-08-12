function els = calc_tetra_childs(els,lc)
%
% Purpose: Calculate child element (1 to 8) inside which the points defined
%          by local coordinates "lc" are located
%
% Input:
%   els    :: parent elements
%   lc     :: local coordinates (r,s,t,u=1-r-s-t) of points in parent elements
% Output:
%   els    :: child elements
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH April 2011
%
    
sub1    = lc(1,:) >= 0.5; % child is subelement #1
sub2    = lc(2,:) >= 0.5; % child is subelement #2
sub3    = lc(3,:) >= 0.5; % child is subelement #3
sub4    = sum(lc)  < 0.5; % child is subelement #4
sub5to8 = ~sub1 & ~sub2 & ~sub3 & ~sub4;
sub57   = sub5to8 & sum(lc(1:2,:))>=0.5;
sub58   = sub5to8 & sum(lc(2:3,:))<=0.5;
sub5    = sub57 & sub58; % child is subelement #5
sub7    = sub57 & ~sub5; % child is subelement #7
sub8    = sub58 & ~sub5; % child is subelement #8

subs       = 6*ones(size(els));
subs(sub1) = 1;
subs(sub2) = 2;
subs(sub3) = 3;
subs(sub4) = 4;
subs(sub5) = 5;
subs(sub7) = 7;
subs(sub8) = 8;

if isinteger(els)
    els = uint32(8*(els-1)) + uint32(subs);
else
    els = 8*(els-1) + subs;
end

% % Old version using indices instead of logical vectors:
% nn   = length(els);
% ind1 = find( LCOORD(1,:) >= 0.5 );
% ind2 = find( LCOORD(2,:) >= 0.5 );
% ind3 = find( LCOORD(3,:) >= 0.5 );
% ind4 = find( sum(LCOORD)  < 0.5 );
% 
% lc5to8  = LCOORD; lc5to8(:,[ind1 ind2 ind3 ind4]) = [];
% ind5to8 = setdiff(1:nn,[ind1 ind2 ind3 ind4]);
% ind57   = ind5to8(sum(lc5to8(1:2,:))>=0.5);
% ind58   = ind5to8(sum(lc5to8(2:3,:))<=0.5);
% 
% ind5 = intersect(ind57,ind58);
% ind7 = setdiff  (ind57,ind5);
% ind8 = setdiff  (ind58,ind5);
% ind6 = setdiff  (1:nn,[ind1 ind2 ind3 ind4 ind5 ind7 ind8]);
% 
% inds       = zeros(size(els));
% inds(ind1) = 1;
% inds(ind2) = 2;
% inds(ind3) = 3;
% inds(ind4) = 4;
% inds(ind5) = 5;
% inds(ind6) = 6;
% inds(ind7) = 7;
% inds(ind8) = 8;

% els = 8*(els-1) + subs;    

end % END OF SUBFUNCTION calc_tetra_childs