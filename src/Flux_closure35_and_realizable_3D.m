function [Fx,Fy,Fz,M4r] = Flux_closure35_and_realizable_3D(M4,flag2D,Ma)
% Flux_closure35_and_realizable_3D computes 3-D fluxes for all moments and
% corrects unrealizable moments
%
%   Input:
%       M4 = 35 moments up to 4th order 
%          = [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,...
%             M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,...
%             M031,M012,M112,M013,M022]
%       flagreal = 1 to check realizability
%   Output:
%       Fx = x flux moments found from closure
%          = [M100,M200,M300,M400,M500,M110,M210,M310,M410,M120,M220,M320,M130,M230,M140,...
%             M101,M201,M301,M401,M102,M202,M302,M103,M203,M104,M111,M211,M311,M121,M221,...
%             M131,M112,M212,M113,M122]
%       Fy = y flux moments found from closure
%          = [M010,M110,M210,M310,M410,M020,M120,M220,M320,M030,M130,M230,M040,M140,M050,...
%             M011,M111,M211,M311,M012,M112,M212,M013,M113,M014,M021,M121,M221,M031,M131,...
%             M041,M022,M122,M023,M032]
%       Fz = z flux moments found from closure
%          = [M001,M101,M201,M301,M401,M011,M111,M211,M311,M021,M121,M221,M031,M131,M041,...
%             M002,M102,M202,M302,M003,M103,M203,M004,M104,M005,M012,M112,M212,M022,M122,...
%             M032,M013,M113,M014,M023]
%       M4r= realizable moments

s3max = 4.0 + abs(Ma)/2.0;
h2min = 1d-8;
itrealmax = 6;

%% compute mean velocities
M000 = M4(1);
M100 = M4(2);
M010 = M4(6);
M001 = M4(16);
umean = M100/M000;
vmean = M010/M000;
wmean = M001/M000; 
%
%% compute central and standardized moments from 35 known moments
[C4,S4] = M2CS4_35(M4);
%
C200 = max(eps,C4(3));
C020 = max(eps,C4(10));
C002 = max(eps,C4(20));
%
S300=S4(4);
S400=S4(5);
S110=S4(7);
S210=S4(8);
S310=S4(9);
S120=S4(11);
S220=S4(12);
S030=S4(13);
S130=S4(14);
S040=S4(15);
%
S101=S4(17);
S201=S4(18);
S301=S4(19);
S102=S4(21);
S202=S4(22);
S003=S4(23);
S103=S4(24);
S004=S4(25);
S011=S4(26);
S111=S4(27);
S211=S4(28);
S021=S4(29);
S121=S4(30);
%
S031=S4(31);
S012=S4(32);
S112=S4(33);
S013=S4(34);
S022=S4(35);
%
%% check univariate moments
H200 = S400 - S300^2 - 1;
H020 = S040 - S030^2 - 1;
H002 = S004 - S003^2 - 1;
% check and correct realizability of S400, S040, S004
if H200 <= h2min
    H200 = h2min;
    S400 = H200 + S300^2 + 1;
end
if H020 <= h2min
    H020 = h2min;
    S040 = H020 + S030^2 + 1;
end
if H002 <= h2min
    H002 = h2min;
    S004 = H002 + S003^2 + 1;
end
% set limits on S300 and S030
if S300 < -s3max
    S300 = -s3max;
    S400 = H200 + S300^2 + 1;
elseif S300 > s3max
    S300 = s3max;
    S400 = H200 + S300^2 + 1;
end
if S030 < -s3max
    S030 = -s3max;
    S040 = H020 + S030^2 + 1;
elseif S030 > s3max
    S030 = s3max;
    S040 = H020 + S030^2 + 1;
end
if S003 < -s3max
    S003 = -s3max;
    S004 = H002 + S003^2 + 1;
elseif S003 > s3max
    S003 = s3max;
    S004 = H002 + S003^2 + 1;
end
%
%% 4-order moments: check maximum bounds on S220, S202, S022
A220 = sqrt((H200+S300^2)*(H020+S030^2));
S220 = realizability_engine('S220', S110, S220, A220);
A202 = sqrt((H200+S300^2)*(H002+S003^2));
S202 = realizability_engine('S220', S101, S202, A202);
A022 = sqrt((H020+S030^2)*(H002+S003^2));
S022 = realizability_engine('S220', S011, S022, A022);
%
%% check and correct realizability of 2D moments
[S300r1,S400r1,S110r,S210r,S310r,S120r,S220r,S030r1,S130r,S040r1] = ...
        check2D(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040);
[S300r2,S400r2,S101r,S201r,S301r,S102r,S202r,S003r2,S103r,S004r2] = ...
        check2D(S300,S400,S101,S201,S301,S102,S202,S003,S103,S004);
[S030r3,S040r3,S011r,S021r,S031r,S012r,S022r,S003r3,S013r,S004r3] = ...
    check2D(S030,S040,S011,S021,S031,S012,S022,S003,S013,S004);
%
% store original values
S111r = S111;
S211r = S211;
S121r = S121;
S112r = S112;
S300r = S300;
S030r = S030;
S003r = S003;
S400r = S400;
S040r = S040;
S004r = S004;
%%
% check for non-realizable 2D moments outside the box
R110 = 1-S110^2;
R101 = 1-S101^2;
R011 = 1-S011^2;
%
if R110 <= 0 || R101 <= 0 || R011 <= 0
    % treat cases where one or more 2D 2nd-order moments is non-realizable (corners and edges)
    disp('treat edges')
    edge_corner_correction
    % recheck 2D cross moments
    [S210,S120,S310,S220,S130] = realizable_2D(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040);
    [S201,S301,S102,S202,S103] = realizable_2D(S300,S400,S101,S201,S301,S102,S202,S003,S103,S004);
    [S021,S031,S012,S022,S013] = realizable_2D(S030,S040,S011,S021,S031,S012,S022,S003,S013,S004);
    %
else
    % treat cases with faces or interior of 2nd-order moment space
    [S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
     S101,S201,S301,S102,S202,S003,S103,S004,S011,S111,...
     S211,S021,S121,S031,S012,S112,S013,S022,flag220] = ...
    realizable_3D(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
                  S101,S201,S301,S102,S202,S003,S103,S004,S011,S111,...
                  S211,S021,S121,S031,S012,S112,S013,S022);
    itreal = 0;
    while flag220 == 1 && itreal < itrealmax
        itreal = itreal + 1;
        % repeat if S220 has changed
        [S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
         S101,S201,S301,S102,S202,S003,S103,S004,S011,S111,...
         S211,S021,S121,S031,S012,S112,S013,S022,~] = ...
        realizable_3D(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
                      S101,S201,S301,S102,S202,S003,S103,S004,S011,S111,...
                      S211,S021,S121,S031,S012,S112,S013,S022);
    end
end
%% force moments for pure 2-D case with S011=S101=0
if flag2D == 1
    if abs(S011) > h2min || abs(S101) > h2min
        warning('flag2D==1 but S011 or S101 is nonzero!')
    end
    S011=0;
    S101=0;
    S201=0;
    S021=0;
    S012=0;
    S102=0;
    S111=0;
    %
    S202=1;
    S022=1;
    S112=S110;
    S121=0;
    S211=0;
    S103=0;
    S013=0;
    S301=0;
    S031=0;
end
% 
%% 3-D hyqmom closures for 5th-order standardized moments
[S500,S410,S320,S230,S140,S401,S302,S203,S104,S311,...
 S221,S131,S212,S113,S122,S050,S041,S032,S023,S014,S005] = ...
hyqmom_3D(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
          S101,S201,S301,S102,S202,S003,S103,S004,S011,S111,...
          S211,S021,S121,S031,S012,S112,S013,S022);
%
%% 5th-order central moments from corrected standardized moments
sC200 = sqrt(max(eps,C200));
sC020 = sqrt(max(eps,C020));
sC002 = sqrt(max(eps,C002));
%
C110=S110*sC200*sC020;
C101=S101*sC200*sC002;
C011=S011*sC020*sC002;
%
C300=S300*sC200^3;
C210=S210*sC200^2*sC020;
C201=S201*sC200^2*sC002;
C120=S120*sC200*sC020^2;
C111=S111*sC200*sC020*sC002;
C102=S102*sC200*sC002^2;
C030=S030*sC020^3;
C021=S021*sC020^2*sC002;
C012=S012*sC020*sC002^2;
C003=S003*sC002^3;
%
C400=S400*sC200^4;
C310=S310*sC200^3*sC020;
C301=S301*sC200^3*sC002;
C220=S220*sC200^2*sC020^2;
C211=S211*sC200^2*sC020*sC002;
C202=S202*sC200^2*sC002^2;
C130=S130*sC200*sC020^3;
C121=S121*sC200*sC020^2*sC002;
C112=S112*sC200*sC020*sC002^2;
C103=S103*sC200*sC002^3;
C040=S040*sC020^4;
C031=S031*sC020^3*sC002;
C022=S022*sC020^2*sC002^2;
C013=S013*sC020*sC002^3;
C004=S004*sC002^4;
%
% closure moments
C500 = S500*sC200^5;
C410 = S410*sC200^4*sC020;
C401 = S401*sC200^4*sC002;
C320 = S320*sC200^3*sC020^2;
C311 = S311*sC200^3*sC020*sC002;
C302 = S302*sC200^3*sC002^2;
C230 = S230*sC200^2*sC020^3;
C221 = S221*sC200^2*sC020^2*sC002;
C212 = S212*sC200^2*sC020*sC002^2;
C203 = S203*sC200^2*sC002^3;
C140 = S140*sC200*sC020^4;
C131 = S131*sC200*sC020^3*sC002;
C122 = S122*sC200*sC020^2*sC002^2;
C113 = S113*sC200*sC020*sC002^3;
C104 = S104*sC200*sC002^4;
C050 = S050*sC020^5;
C041 = S041*sC020^4*sC002;
C032 = S032*sC020^3*sC002^2;
C023 = S023*sC020^2*sC002^3;
C014 = S014*sC020*sC002^4;
C005 = S005*sC002^5;
%
%% 5th-order moments from central moments
M5 = C5toM5_3D(M000,umean,vmean,wmean,C200,C110,C101,C020,C011,C002,...
               C300,C210,C201,C120,C111,C102,C030,C021,C012,C003,...
               C400,C310,C301,C220,C211,C202,C130,C121,C112,C103,C040,C031,C022,C013,C004,...
               C500,C410,C320,C230,C140,C401,C302,C203,C104,C311,C221,C131,C212,C113,C122,C050,C041,C032,C023,C014,C005);
%
M000 = M5(1,1,1);
M100 = M5(2,1,1);
M010 = M5(1,2,1);
M001 = M5(1,1,2);
%
M200 = M5(3,1,1);
M110 = M5(2,2,1);
M101 = M5(2,1,2);
M020 = M5(1,3,1);
M011 = M5(1,2,2);
M002 = M5(1,1,3);
%
M300 = M5(4,1,1);
M210 = M5(3,2,1);
M201 = M5(3,1,2);
M120 = M5(2,3,1);
M111 = M5(2,2,2);
M102 = M5(2,1,3);
M030 = M5(1,4,1);
M021 = M5(1,3,2);
M012 = M5(1,2,3);
M003 = M5(1,1,4);
%
M400 = M5(5,1,1);
M310 = M5(4,2,1);
M301 = M5(4,1,2);
M220 = M5(3,3,1);
M211 = M5(3,2,2);
M202 = M5(3,1,3);
M130 = M5(2,4,1);
M121 = M5(2,3,2);
M112 = M5(2,2,3);
M103 = M5(2,1,4);
M040 = M5(1,5,1);
M031 = M5(1,4,2);
M022 = M5(1,3,3);
M013 = M5(1,2,4);
M004 = M5(1,1,5);
%
% closure moments
M500 = M5(6,1,1);
M410 = M5(5,2,1); 
M320 = M5(4,3,1);
M230 = M5(3,4,1);
M140 = M5(2,5,1);
M401 = M5(5,1,2); 
M302 = M5(4,1,3); 
M203 = M5(3,1,4); 
M104 = M5(2,1,5); 
M311 = M5(4,2,2); 
M221 = M5(3,3,2); 
M131 = M5(2,4,2); 
M212 = M5(3,2,3); 
M113 = M5(2,2,4); 
M122 = M5(2,3,3); 
M050 = M5(1,6,1); 
M041 = M5(1,5,2); 
M032 = M5(1,4,3); 
M023 = M5(1,3,4); 
M014 = M5(1,2,5); 
M005 = M5(1,1,6);
%
%% flux closures
Fx = [M100,M200,M300,M400,M500,M110,M210,M310,M410,M120,M220,M320,M130,M230,M140,...
      M101,M201,M301,M401,M102,M202,M302,M103,M203,M104,M111,M211,M311,M121,M221,...
      M131,M112,M212,M113,M122];
Fy = [M010,M110,M210,M310,M410,M020,M120,M220,M320,M030,M130,M230,M040,M140,M050,...
      M011,M111,M211,M311,M012,M112,M212,M013,M113,M014,M021,M121,M221,M031,M131,...
      M041,M022,M122,M023,M032];
Fz = [M001,M101,M201,M301,M401,M011,M111,M211,M311,M021,M121,M221,M031,M131,M041,...
      M002,M102,M202,M302,M003,M103,M203,M004,M104,M005,M012,M112,M212,M022,M122,...
      M032,M013,M113,M014,M023];
%
% realizable moments
M4r= [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,...
      M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,...
      M031,M012,M112,M013,M022];
%
end