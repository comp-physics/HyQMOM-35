function m = unpack_M(M)
% unpack_M Convert moment vector to struct with named fields
%
% Input:
%   M - 35-element moment vector [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,...
%                                 M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,...
%                                 M031,M012,M112,M013,M022]
%
% Output:
%   m - struct with named moment fields

m = struct();

% 0th order
m.m000 = M(1);

% 1st order
m.m100 = M(2);
m.m010 = M(6);
m.m001 = M(16);

% 2nd order
m.m200 = M(3);
m.m110 = M(7);
m.m101 = M(17);
m.m020 = M(10);
m.m011 = M(26);
m.m002 = M(20);

% 3rd order
m.m300 = M(4);
m.m210 = M(8);
m.m201 = M(18);
m.m120 = M(11);
m.m111 = M(27);
m.m102 = M(21);
m.m030 = M(13);
m.m021 = M(29);
m.m012 = M(32);
m.m003 = M(23);

% 4th order
m.m400 = M(5);
m.m310 = M(9);
m.m301 = M(19);
m.m220 = M(12);
m.m211 = M(28);
m.m202 = M(22);
m.m130 = M(14);
m.m121 = M(30);
m.m112 = M(33);
m.m103 = M(24);
m.m040 = M(15);
m.m031 = M(31);
m.m022 = M(35);
m.m013 = M(34);
m.m004 = M(25);

end
