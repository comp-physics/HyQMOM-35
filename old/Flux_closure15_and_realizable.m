function [Fx,Fy,M4r] = Flux_closure15_and_realizable(M4,flagreal,Ma)
% Flux_closure15_and_realizable computes 2-D fluxes for all moments and
% corrects moments
%
%   Input:
%       M4 = moments up to 4th order 
%          = [M00,M10,M20,M30,M40,M01,M11,M21,M31,M02,M12,M22,M03,M13,M04]
%       flagreal = 1 to check realizability
%   Output:
%       Fx = x flux moments found from closure
%          = [M10 M20 M30 M40 M50 M11 M21 M31 M41 M12 M22 M32 M13 M23 M14]
%       Fy = y flux moments found from closure
%          = [M01 M11 M21 M31 M41 M02 M12 M22 M32 M03 M13 M23 M04 M14 M05]
%       M4r= realizable moments

s3max = 4 + abs(Ma)/2;
h2min = 1000*eps;

% compute mean velocities
M00 = M4(1);
M10 = M4(2);
M01 = M4(6);
umean = M10/M00;
vmean = M01/M00;
%
% compute central and standardized moments and closures from 15 known moments
[C4,S4] = M2CS4_15(M4);
%
C20 = C4(3);
C02 = C4(10);
%
S11 = S4(7);
S30 = S4(4);
S21 = S4(8);
S12 = S4(11);
S03 = S4(13);
S40 = S4(5);
S31 = S4(9);
S22 = S4(12);
S13 = S4(14);
S04 = S4(15);
%
H20 = S40 - S30^2 - 1;
H02 = S04 - S03^2 - 1;
% check and correct realizability of S40 and S04
if H20 <= h2min
    H20 = h2min;
    S40 = H20 + S30^2 + 1;
end
if H02 <= h2min
    H02 = h2min;
    S04 = H02 + S03^2 + 1;
end
% set limits on S30 and S03
if S30 < -s3max
    S30 = -s3max;
    S40 = H20 + S30^2 + 1;
elseif S30 > s3max
    S30 = s3max;
    S40 = H20 + S30^2 + 1;
end
if S03 < -s3max
    S03 = -s3max;
    S04 = H02 + S03^2 + 1;
elseif S03 > s3max
    S03 = s3max;
    S04 = H02 + S03^2 + 1;
end
%
% check and correct realizability of S11
if S11 >= 1
    S11 = 1;
    s3m = sqrt(abs(S30*S03));
    s4m = sqrt(S40*S04);
    H2m = max(h2min,s4m - s3m^2 - 1);
    s4m = H2m + s3m^2 + 1;
    S30 = sign(S30)*s3m;
    S03 = S11*S30;
    S40 = s4m;
    S04 = s4m;
    S12 = S11*S03;
    S21 = S11*S30;
    S13 = S11*S04;
    S31 = S11*S40;
    S22 = s4m;
elseif S11 <= -1
    S11 = -1;
    s3m = sqrt(abs(S30*S03));
    s4m = sqrt(S40*S04);
    H2m = max(h2min,s4m - s3m^2 - 1);
    s4m = H2m + s3m^2 + 1;
    S30 = sign(S30)*s3m;
    S03 = S11*S30;
    S40 = s4m;
    S04 = s4m;
    S12 = S11*S03;
    S21 = S11*S30;
    S13 = S11*S04;
    S31 = S11*S40;
    S22 = s4m;
else
    %
    % check and correct realizability of cross moments
    if flagreal == 1
        [S21,S12,S31,S22,S13] = realizable_2D(S30,S40,S11,S21,S31,S12,S22,S03,S13,S04);
    end
end
% hyqmom closures for 5th order standardized moments
[S50,S41,S32,S23,S14,S05] = hyqmom_2D(S30,S40,S11,S21,S31,S12,S22,S03,S13,S04);
%
% central moments from standized moments
sC20= sqrt(C20);
sC02= sqrt(C02);
C11 = S11*sC20*sC02;
C30 = S30*sC20^3; 
C21 = S21*sC20^2*sC02;
C12 = S12*sC20*sC02^2;
C03 = S03*sC02^3;
C40 = S40*sC20^4; 
C31 = S31*sC20^3*sC02; 
C22 = S22*sC20^2*sC02^2;
C13 = S13*sC20*sC02^3;
C04 = S04*sC02^4;
C50 = S50*sC20^5;
C41 = S41*sC20^4*sC02; 
C32 = S32*sC20^3*sC02^2; 
C23 = S23*sC20^2*sC02^3;
C14 = S14*sC20*sC02^4;
C05 = S05*sC02^5;
%
% integer moments from central moments
M5 = C5toM5(M00,umean,vmean,C20,C11,C02,C30,C21,C12,C03,C40,C31,C22,C13,C04,C50,C41,C32,C23,C14,C05);
%
M00 = M5(1,1);
M10 = M5(2,1);
M01 = M5(1,2);
M20 = M5(3,1);
M11 = M5(2,2);
M02 = M5(1,3);
M30 = M5(4,1);
M21 = M5(3,2);
M12 = M5(2,3);
M03 = M5(1,4);
M40 = M5(5,1); 
M31 = M5(4,2);
M22 = M5(3,3);
M13 = M5(2,4);
M04 = M5(1,5);
M50 = M5(6,1);
M41 = M5(5,2); 
M32 = M5(4,3);
M23 = M5(3,4);
M14 = M5(2,5);
M05 = M5(1,6);
%
% flux closures
Fx = [M10 M20 M30 M40 M50 M11 M21 M31 M41 M12 M22 M32 M13 M23 M14];
Fy = [M01 M11 M21 M31 M41 M02 M12 M22 M32 M03 M13 M23 M04 M14 M05];
%
% realizable moments
M4r= [M00,M10,M20,M30,M40,M01,M11,M21,M31,M02,M12,M22,M03,M13,M04];
end