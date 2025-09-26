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

% Get centralized constants
constants = get_constants(Ma);
s3max = constants.s3max;
h2min = constants.h2min;
itrealmax = constants.itrealmax;

%% compute mean velocities
[umean, vmean, wmean] = means_from_M(M4);
idx = moment_indices();
M000 = M4(idx.m000); 
%
%% compute central and standardized moments from 35 known moments
[C4,S4] = M2CS4_35(M4);
%
C200 = max(eps,C4(idx.C200));
C020 = max(eps,C4(idx.C020));
C002 = max(eps,C4(idx.C002));
%
% Extract standardized moments using indices
S300=S4(idx.S300); S400=S4(idx.S400); S030=S4(idx.S030); S040=S4(idx.S040); S003=S4(idx.S003); S004=S4(idx.S004);
S110=S4(7); S210=S4(8); S310=S4(9); S120=S4(11); S220=S4(12); S130=S4(14);
S101=S4(17); S201=S4(18); S301=S4(19); S102=S4(21); S202=S4(22); S103=S4(24);
S011=S4(26); S111=S4(27); S211=S4(28); S021=S4(29); S121=S4(30);
S031=S4(31); S012=S4(32); S112=S4(33); S013=S4(34); S022=S4(35);
%
%% check univariate moments
[S300, S400, H200] = enforce_univariate(S300, S400, h2min, s3max);
[S030, S040, H020] = enforce_univariate(S030, S040, h2min, s3max);
[S003, S004, H002] = enforce_univariate(S003, S004, h2min, s3max);
%
%% 4-order moments: check maximum bounds on S220, S202, S022
S220 = cap_S220(S110, S220, H200, H020, S300, S030);
S202 = cap_S220(S101, S202, H200, H002, S300, S003);
S022 = cap_S220(S011, S022, H020, H002, S030, S003);
%
%% check and correct realizability of 2D moments
% Apply check2D to all three coordinate planes
S4_checked = check2D_all_planes(S4);

% Update local variables with checked values
S110r = S4_checked(7); S210r = S4_checked(8); S310r = S4_checked(9); 
S120r = S4_checked(11); S220r = S4_checked(12); S130r = S4_checked(14);
S101r = S4_checked(17); S201r = S4_checked(18); S301r = S4_checked(19);
S102r = S4_checked(21); S202r = S4_checked(22); S103r = S4_checked(24);
S011r = S4_checked(26); S021r = S4_checked(29); S031r = S4_checked(31);
S012r = S4_checked(32); S013r = S4_checked(34); S022r = S4_checked(35);

% Store original values of moments that aren't corrected by check2D
S111r = S111; S211r = S211; S121r = S121; S112r = S112;
S300r = S4_checked(4); S030r = S4_checked(13); S003r = S4_checked(23);
S400r = S4_checked(5); S040r = S4_checked(15); S004r = S4_checked(25);
%%
% check for non-realizable 2D moments outside the box
R110 = 1-S110^2;
R101 = 1-S101^2;
R011 = 1-S011^2;
%
if R110 <= 0 || R101 <= 0 || R011 <= 0
    % treat cases where one or more 2D 2nd-order moments is non-realizable (corners and edges)
    fprintf('Treating edge/corner case: R110=%g, R101=%g, R011=%g\n', R110, R101, R011);
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
[sC200, sC020, sC002] = sC_from_C4(C4);
% Scale standardized moments to central moments
% 2nd order
C110=S110*sC200*sC020; C101=S101*sC200*sC002; C011=S011*sC020*sC002;
% 3rd order  
C300=S300*sC200^3; C210=S210*sC200^2*sC020; C201=S201*sC200^2*sC002;
C120=S120*sC200*sC020^2; C111=S111*sC200*sC020*sC002; C102=S102*sC200*sC002^2;
C030=S030*sC020^3; C021=S021*sC020^2*sC002; C012=S012*sC020*sC002^2; C003=S003*sC002^3;
% 4th order
C400=S400*sC200^4; C310=S310*sC200^3*sC020; C301=S301*sC200^3*sC002;
C220=S220*sC200^2*sC020^2; C211=S211*sC200^2*sC020*sC002; C202=S202*sC200^2*sC002^2;
C130=S130*sC200*sC020^3; C121=S121*sC200*sC020^2*sC002; C112=S112*sC200*sC020*sC002^2; C103=S103*sC200*sC002^3;
C040=S040*sC020^4; C031=S031*sC020^3*sC002; C022=S022*sC020^2*sC002^2; C013=S013*sC020*sC002^3; C004=S004*sC002^4;

% 5th order closure moments using helper
S5_cell = {S500,S410,S320,S230,S140,S401,S302,S203,S104,S311,S221,S131,S212,S113,S122,S050,S041,S032,S023,S014,S005};
C5_cell = scale_S5_to_C5(S5_cell, sC200, sC020, sC002);
[C500,C410,C320,C230,C140,C401,C302,C203,C104,C311,C221,C131,C212,C113,C122,C050,C041,C032,C023,C014,C005] = C5_cell{:};
%
%% 5th-order moments from central moments
M5 = C5toM5_3D(M000,umean,vmean,wmean,C200,C110,C101,C020,C011,C002,...
               C300,C210,C201,C120,C111,C102,C030,C021,C012,C003,...
               C400,C310,C301,C220,C211,C202,C130,C121,C112,C103,C040,C031,C022,C013,C004,...
               C500,C410,C320,C230,C140,C401,C302,C203,C104,C311,C221,C131,C212,C113,C122,C050,C041,C032,C023,C014,C005);
% Extract moments from M5 array
% 0th-1st order
M000 = M5(1,1,1); M100 = M5(2,1,1); M010 = M5(1,2,1); M001 = M5(1,1,2);
% 2nd order
M200 = M5(3,1,1); M110 = M5(2,2,1); M101 = M5(2,1,2); M020 = M5(1,3,1); M011 = M5(1,2,2); M002 = M5(1,1,3);
% 3rd order
M300 = M5(4,1,1); M210 = M5(3,2,1); M201 = M5(3,1,2); M120 = M5(2,3,1); M111 = M5(2,2,2);
M102 = M5(2,1,3); M030 = M5(1,4,1); M021 = M5(1,3,2); M012 = M5(1,2,3); M003 = M5(1,1,4);
% 4th order
M400 = M5(5,1,1); M310 = M5(4,2,1); M301 = M5(4,1,2); M220 = M5(3,3,1); M211 = M5(3,2,2);
M202 = M5(3,1,3); M130 = M5(2,4,1); M121 = M5(2,3,2); M112 = M5(2,2,3); M103 = M5(2,1,4);
M040 = M5(1,5,1); M031 = M5(1,4,2); M022 = M5(1,3,3); M013 = M5(1,2,4); M004 = M5(1,1,5);
% 5th order closure moments
M500 = M5(6,1,1); M410 = M5(5,2,1); M320 = M5(4,3,1); M230 = M5(3,4,1); M140 = M5(2,5,1);
M401 = M5(5,1,2); M302 = M5(4,1,3); M203 = M5(3,1,4); M104 = M5(2,1,5);
M311 = M5(4,2,2); M221 = M5(3,3,2); M131 = M5(2,4,2); M212 = M5(3,2,3); M113 = M5(2,2,4); M122 = M5(2,3,3);
M050 = M5(1,6,1); M041 = M5(1,5,2); M032 = M5(1,4,3); M023 = M5(1,3,4); M014 = M5(1,2,5); M005 = M5(1,1,6);
%
%% flux closures - generate programmatically using indices
% Extract flux moments using subscript indices
Fx = zeros(35,1); Fy = zeros(35,1); Fz = zeros(35,1);
for k = 1:35
    Fx(k) = M5(idx.Fx_subs(k,1), idx.Fx_subs(k,2), idx.Fx_subs(k,3));
    Fy(k) = M5(idx.Fy_subs(k,1), idx.Fy_subs(k,2), idx.Fy_subs(k,3));
    Fz(k) = M5(idx.Fz_subs(k,1), idx.Fz_subs(k,2), idx.Fz_subs(k,3));
end
Fx = Fx'; Fy = Fy'; Fz = Fz';
%
% realizable moments
M4r= [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,...
      M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,...
      M031,M012,M112,M013,M022];
%
end