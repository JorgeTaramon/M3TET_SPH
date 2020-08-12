function [ipx,ipw] = ip_tetrahedron(nip)
% Usage: [ipx,ipw] = ip_tetrahedron(nip)
% 
% Purpose: Defines integration points and weights
%
% Input:
%   nip   : [scalar]    : number of integration points
%
% Output:
%   ipx   : [matrix]    : local coordinates of IPs (3 x nip)
%   ipw   : [rowvector] : integration weights at each IP
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Dec 2012
% JH Mar 2013
% JH Dec 2014: changed orientation of ipx; now [3 x nip]
% JH Mar 2016: started adding higher order integration rules
% JMT Jul 2017: fixed integration rule for 14 points
%

switch nip
    case 1
        % ONE POINT INTEGRATION (degree of precision: 1)
        % Taken from Zienkiewicz Vol 1 (4th) edition p.177
        % (also in Hughes p.174)
        ipx    = zeros(3,nip); % one point defined by the 3 local coords r, s, t
        ipx(1) = 1/4; % r
        ipx(2) = 1/4; % s
        ipx(3) = 1/4; % t
        ipw    = 1;   % weight for integration
        
    case 4
        % FOUR POINT INTEGRATION (degree of precision: 2)
        % Taken from Zienkiewicz Vol 1 (4th) edition p.177 (also in Hughes p.174)
        a        = (5+3*sqrt(5))/20; % 0.58541020;
        b        = (5-sqrt(5))/20;   % 0.13819660;
        ipx      = zeros(3,nip); % 4 points defined by the 3 local coords r, s, t
        ipx(1,:) = [ a   b   b   b ]; % r at all points
        ipx(2,:) = [ b   a   b   b ]; % s at all points
        ipx(3,:) = [ b   b   a   b ]; % t at all points
        ipw      = [1/4 1/4 1/4 1/4]; % weights for integration

    case 5
        % FIVE POINT INTEGRATION (degree of precision: 3)
        % Taken from Zienkiewicz Vol 1 (4th) edition p.177 (also in Hughes p.174)
        ipx      = zeros(3,nip); % 5 points defined by the 3 local coords r, s, t
        ipx(1,:) = [ 1/4  1/2  1/6  1/6  1/6]; % r at all points
        ipx(2,:) = [ 1/4  1/6  1/2  1/6  1/6]; % s at all points
        ipx(3,:) = [ 1/4  1/6  1/6  1/2  1/6]; % t at all points
        ipw      = [-4/5 9/20 9/20 9/20 9/20]; % weights for integration
        
    case 8
        % EIGHT POINT INTEGRATION (degree of precision: 3)
        % Taken from Table III in Gellert & Harbord, (1991), Moderate degree
        % cubature formulas for 3-D tetrahedral finite-element approximations.
        % Communications in Applied Numerical Methods, Vol.7, 487-495.
        ipx      = zeros(3,nip); % 8 points defined by the 3 local coords r, s, t
        ipx(1,:) = [  0   1/3  1/3  1/3   1    0    0    0 ]; % r at all points
        ipx(2,:) = [ 1/3   0   1/3  1/3   0    1    0    0 ]; % s at all points
        ipx(3,:) = [ 1/3  1/3   0   1/3   0    0    1    0 ]; % t at all points
        ipw      = [9/40 9/40 9/40 9/40 1/40 1/40 1/40 1/40]; % weights for integration

    case 10
        % TEN POINT INTEGRATION (degree of precision: 3)
        % Taken from Appendix F in 
        a  = 0.7784952948213300;
        b  = 0.0738349017262234;
        c  = 0.4062443438840510;
        d  = 0.0937556561159491;
        w1 = 0.0476331348432089;
        w2 = 0.1349112434378610;
        ipx      = zeros(3,nip); % 10 points defined by the 3 local coords r, s, t
        ipx(1,:) = [ a  b  b  b  c  c  c  d  d  d]; % r at all points
        ipx(2,:) = [ b  a  b  b  c  d  d  c  c  d]; % s at all points
        ipx(3,:) = [ b  b  a  b  d  c  d  c  d  c]; % t at all points
        ipw      = [w1 w1 w1 w1 w2 w2 w2 w2 w2 w2]; % weights for integration
        
    case 11
        % ELEVEN POINT (degree of precision: 4)
        % Taken from Table IV in Gellert & Harbord, (1991), Moderate degree
        % cubature formulas for 3-D tetrahedral finite-element approximations.
        % Communications in Applied Numerical Methods, Vol.7, 487-495.
        a        = (1 + sqrt(5/14))/4;
        b        = (1 - sqrt(5/14))/4;
        w1       = -148/1875;
        w2       =  343/7500;
        w3       =   56/375;
        ipx      = zeros(3,nip); % 11 points defined by the 3 local coords r, s, t
        ipx(1,:) = [1/4 11/14  1/14  1/14 1/14  a  a  a  b  b  b]; % r at all points
        ipx(2,:) = [1/4  1/14 11/14  1/14 1/14  a  b  b  a  a  b]; % s at all points
        ipx(3,:) = [1/4  1/14  1/14 11/14 1/14  b  a  b  a  b  a]; % t at all points
        ipw      = [ w1   w2    w2    w2   w2  w3 w3 w3 w3 w3 w3]; % weights for integration
        
    case 14
        % FOURTEEN POINT (degree of precision: 5)
        % Taken from Table VI in Gellert & Harbord, (1991), Moderate degree
        % cubature formulas for 3-D tetrahedral finite-element approximations.
        % Communications in Applied Numerical Methods, Vol.7, 487-495.
        g = 1/(46*sqrt(46));
        h = acos(g) + (2/3)*asin(g);
        k = (104 + 8*sqrt(46)*cos(h))/3;
        s = sqrt(49 - k);
        b = (7 + s)/k; % = 0.3108859192633005
        a = 1 - 3*b;   % = 0.0673422422100983
        d = (7 - s)/k; % = 0.0927352503108912
        c = 1 - 3*d;   % = 0.7217942490673264
        p = (98 - k - 14*s)/(1680*s*(b - a)^3); % = 0.1126879257180162
        q = (98 - k + 14*s)/(1680*s*(c - d)^3); % = 0.0734930431163619
        r = (1 - 4*(p + q))/6;                  % = 0.0425460207770812
        e = (1 + (2/(105*r))^(1/4))/4; % = 0.4544962958743506
        f = (1 - 2*e)/2;               % = 0.0455037041256494
        ipx      = zeros(3,nip); % 14 points defined by the 3 local coords r, s, t
        ipx(1,:) = [a b b b c d d d e e e f f f]; % r at all points
        ipx(2,:) = [b a b b d c d d e f f e e f]; % s at all points
        ipx(3,:) = [b b a b d d c d f e f e f e]; % t at all points
        ipw      = [p p p p q q q q r r r r r r]; % weights for integration

    case 20
        error('This integration rule does not work (b,a,d,etc are not the right values)');
        % TWENTY POINT INTEGRATION (degree of precision: 5)
        % Taken from Zienkiewicz Vol 1 (4th) edition p.177 (also in Hughes p.174)
        % Taken from Table VII in Gellert & Harbord, (1991), Moderate degree
        % cubature formulas for 3-D tetrahedral finite-element approximations.
        % Communications in Applied Numerical Methods, Vol.7, 487-495.
        w1 = 81/2240;
        w2 = 161051/2304960;
        w3 = 409/31395;
        w4 = 2679769/32305455;
        a  = (1 + 21/sqrt(1637))/4; % = 0.3797582452067875
        b  = (1 - 21/sqrt(1637))/4; % = 0.1202417547932126
        ipx      = zeros(3,nip); % 5 points defined by the 3 local coords r, s, t
        ipx(1,:) = [ 0  1/3 1/3 1/3 8/11 1/11 1/11 1/11  0   0   0  1/2 1/2 1/2  a  a  a  b  b  b]; % r at all points
        ipx(2,:) = [1/3  0  1/3 1/3 1/11 8/11 1/11 1/11  0  1/2 1/2  0  1/2 1/2  a  b  b  a  a  b]; % s at all points
        ipx(3,:) = [1/3 1/3  0  1/3 1/11 1/11 8/11 1/11 1/2  0  1/2  0  1/2  0   b  a  b  a  b  a]; % t at all points
        ipw      = [ w1  w1  w1  w1  w2   w2   w2   w2   w3  w3  w3  w3  w3  w3 w4 w4 w4 w4 w4 w4]; % weights for integration
%         % Table VII as in the paper: ONLY 19 points !!!
%         ipx(1,:) = [ 0  1/3 1/3 1/3 8/11 1/11 1/11 1/11  0   0   0  1/2 1/2  a  a  a  b  b  b]; % r at all points
%         ipx(2,:) = [1/3  0  1/3 1/3 1/11 8/11 1/11 1/11  0  1/2 1/2  0  1/2  a  b  b  a  a  b]; % s at all points
%         ipx(3,:) = [1/3 1/3  0  1/3 1/11 1/11 8/11 1/11 1/2  0  1/2  0   0   b  a  b  a  b  a]; % t at all points
%         ipw      = [ w1  w1  w1  w1  w2   w2   w2   w2   w3  w3  w3  w3  w3 w4 w4 w4 w4 w4 w4]; % weights for integration
        
    case 35
        error('This integration rule is not yet coded.');
% 1 0.9197896733368800 0.0267367755543735 0.0267367755543735 0.0267367755543735 0.0021900463965388
% 2 0.0267367755543735 0.9197896733368800 0.0267367755543735 0.0267367755543735 0.0021900463965388
% 3 0.0267367755543735 0.0267367755543735 0.9197896733368800 0.0267367755543735 0.0021900463965388
% 4 0.0267367755543735 0.0267367755543735 0.0267367755543735 0.9197896733368800 0.0021900463965388
% 5 0.1740356302468940 0.7477598884818090 0.0391022406356488 0.0391022406356488 0.0143395670177665
% 6 0.7477598884818090 0.1740356302468940 0.0391022406356488 0.0391022406356488 0.0143395670177665
% 7 0.1740356302468940 0.0391022406356488 0.7477598884818090 0.0391022406356488 0.0143395670177665
% 8 0.7477598884818090 0.0391022406356488 0.1740356302468940 0.0391022406356488 0.0143395670177665
% 9 0.1740356302468940 0.0391022406356488 0.0391022406356488 0.7477598884818090 0.0143395670177665
% 10 0.7477598884818090 0.0391022406356488 0.0391022406356488 0.1740356302468940 0.0143395670177665
% 11 0.0391022406356488 0.1740356302468940 0.7477598884818090 0.0391022406356488 0.0143395670177665
% 12 0.0391022406356488 0.7477598884818090 0.1740356302468940 0.0391022406356488 0.0143395670177665
% 13 0.0391022406356488 0.1740356302468940 0.0391022406356488 0.7477598884818090 0.0143395670177665
% 14 0.0391022406356488 0.7477598884818090 0.0391022406356488 0.1740356302468940 0.0143395670177665
% 15 0.0391022406356488 0.0391022406356488 0.1740356302468940 0.7477598884818090 0.0143395670177665
% 16 0.0391022406356488 0.0391022406356488 0.7477598884818090 0.1740356302468940 0.0143395670177665
% 17 0.4547545999844830 0.4547545999844830 0.0452454000155172 0.0452454000155172 0.0250305395686746
% 18 0.4547545999844830 0.0452454000155172 0.4547545999844830 0.0452454000155172 0.0250305395686746
% 19 0.4547545999844830 0.0452454000155172 0.0452454000155172 0.4547545999844830 0.0250305395686746
% 20 0.0452454000155172 0.4547545999844830 0.4547545999844830 0.0452454000155172 0.0250305395686746
% 21 0.0452454000155172 0.4547545999844830 0.0452454000155172 0.4547545999844830 0.0250305395686746
% 22 0.0452454000155172 0.0452454000155172 0.4547545999844830 0.4547545999844830 0.0250305395686746
% 23 0.5031186450145980 0.2232010379623150 0.2232010379623150 0.0504792790607720 0.0479839333057554
% 24 0.2232010379623150 0.5031186450145980 0.2232010379623150 0.0504792790607720 0.0479839333057554
% 25 0.2232010379623150 0.2232010379623150 0.5031186450145980 0.0504792790607720 0.0479839333057554
% 26 0.5031186450145980 0.2232010379623150 0.0504792790607720 0.2232010379623150 0.0479839333057554
% 27 0.2232010379623150 0.5031186450145980 0.0504792790607720 0.2232010379623150 0.0479839333057554
% 28 0.2232010379623150 0.2232010379623150 0.0504792790607720 0.5031186450145980 0.0479839333057554
% 29 0.5031186450145980 0.0504792790607720 0.2232010379623150 0.2232010379623150 0.0479839333057554
% 30 0.2232010379623150 0.0504792790607720 0.5031186450145980 0.2232010379623150 0.0479839333057554
% 31 0.2232010379623150 0.0504792790607720 0.2232010379623150 0.5031186450145980 0.0479839333057554
% 32 0.0504792790607720 0.5031186450145980 0.2232010379623150 0.2232010379623150 0.0479839333057554
% 33 0.0504792790607720 0.2232010379623150 0.5031186450145980 0.2232010379623150 0.0479839333057554
% 34 0.0504792790607720 0.2232010379623150 0.2232010379623150 0.5031186450145980 0.0479839333057554
% 35 0.2500000000000000 0.2500000000000000 0.2500000000000000 0.2500000000000000 0.0931745731195340
    otherwise
        error('Unkown integration rule.');
end

% Volume of tetrahedron is 1/6 of encompassing cube
% V=0.5*A*h, A=area of base, h=height
ipw = ipw ./ 6;

end % END OF FUNCTION ip_tetrahedron