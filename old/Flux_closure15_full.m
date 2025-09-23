function [Fx,Fy] = Flux_closure15_full(M4)
% Flux_closure15_full computes 2-D fluxes for all moments
%
%   Input:
%       M4 = moments up to 4th order 
%          = [M00,M10,M20,M30,M40,M01,M11,M21,M31,M02,M12,M22,M03,M13,M04]
%   Output:
%       Fx = x flux moments found from closure
%          = [M10 M20 M30 M40 M50 M11 M21 M31 M41 M12 M22 M32 M13 M23 M14]
%       Fy = y flux moments found from closure
%          = [M01 M11 M21 M31 M41 M02 M12 M22 M32 M03 M13 M23 M04 M14 M05]

% compute mean velocities
M00 = M4(1);
M10 = M4(2);
M01 = M4(6);
umean = M10/M00;
vmean = M01/M00;

% compute central and standardized moments and closures from 10 known moments
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
% check realizability of S11
if S11^2 > 1
    S11 = sign(S11);
end
% check realizability of S12 and S21
Del1= max(0,1 - S11^2);
H20 = max(0,S40 - S30^2 - 1);
H02 = max(0,S04 - S03^2 - 1);
Low = S11*S03 - sqrt(Del1*H02);
High= S11*S03 + sqrt(Del1*H02);
if S12 < Low
    S12=Low;
elseif S12 > High
    S12=High;
end
Low = S11*S30 - sqrt(Del1*H20);
High= S11*S30 + sqrt(Del1*H20);
if S21 < Low
    S21=Low;
elseif S21 > High
    S21=High;
end
%
% hyqmom closures
S50 = 0.5*S30*(5*S40 - 3*S30^2 - 1);
S05 = 0.5*S03*(5*S04 - 3*S03^2 - 1);
S32 = ((S03*S30^2)/2 - S03 + S03*(S30^2 - S40 + 1))*S11 + (S40 - (3*S30^2)/2)*S12 + (-(3*S03*S30)/2)*S21 + ((3*S30)/2)*S22 + S03*S31 - S30/2;
S23 = ((S30*S03^2)/2 - S30 + S30*(S03^2 - S04 + 1))*S11 + (S04 - (3*S03^2)/2)*S21 + (-(3*S30*S03)/2)*S12 + ((3*S03)/2)*S22 + S30*S13 - S03/2;
S41 = (S30 - 2*S30*S40 + (9*S30^3)/4)*S11 + ((5*S40)/2 - (15*S30^2)/4 - 3/2)*S21 + (2*S30)*S31;
S14 = (S03 - 2*S03*S04 + (9*S03^3)/4)*S11 + ((5*S04)/2 - (15*S03^2)/4 - 3/2)*S12 + (2*S03)*S13;
%
% closures for central moments from standized moments
sC20 = sqrt(C20);
sC02 = sqrt(C02);
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
% closure moments from central moments
M5 = C5toM5(M00,umean,vmean,C20,C11,C02,C30,C21,C12,C03,C40,C31,C22,C13,C04,C50,C41,C32,C23,C14,C05);
%
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
end