function M = pack_M(m)
% pack_M Convert struct with named moment fields to moment vector
%
% Input:
%   m - struct with named moment fields (m000, m100, etc.)
%
% Output:
%   M - 35-element moment vector [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,...
%                                 M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,...
%                                 M031,M012,M112,M013,M022]

M = zeros(35, 1);

% 0th order
M(1) = m.m000;

% 1st order
M(2) = m.m100;
M(6) = m.m010;
M(16) = m.m001;

% 2nd order
M(3) = m.m200;
M(7) = m.m110;
M(17) = m.m101;
M(10) = m.m020;
M(26) = m.m011;
M(20) = m.m002;

% 3rd order
M(4) = m.m300;
M(8) = m.m210;
M(18) = m.m201;
M(11) = m.m120;
M(27) = m.m111;
M(21) = m.m102;
M(13) = m.m030;
M(29) = m.m021;
M(32) = m.m012;
M(23) = m.m003;

% 4th order
M(5) = m.m400;
M(9) = m.m310;
M(19) = m.m301;
M(12) = m.m220;
M(28) = m.m211;
M(22) = m.m202;
M(14) = m.m130;
M(30) = m.m121;
M(33) = m.m112;
M(24) = m.m103;
M(15) = m.m040;
M(31) = m.m031;
M(35) = m.m022;
M(34) = m.m013;
M(25) = m.m004;

end
