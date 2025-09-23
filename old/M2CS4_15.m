function [C4,S4] = M2CS4_15(M4)
% M(1:Nmom+1): moments of order 0 to Nmom-1
% 
% central and standardized moments
C4 = 0*M4;
S4 = 0*C4;

% rearrange moment vector
M00 = M4(1);
M10 = M4(2);
M01 = M4(6);
M20 = M4(3);
M11 = M4(7);
M02 = M4(10);
M30 = M4(4);
M21 = M4(8);
M12 = M4(11);
M03 = M4(13);
M40 = M4(5);
M31 = M4(9);
M22 = M4(12);
M13 = M4(14);
M04 = M4(15);

% compute central moments
C = M4toC4(M00,M10,M01,M20,M11,M02,M30,M21,M12,M03,M40,M31,M22,M13,M04);

C00 = C(1,1);
C10 = C(2,1);
C01 = C(1,2);
C20 = C(3,1);
C11 = C(2,2);
C02 = C(1,3);
C30 = C(4,1);
C21 = C(3,2);
C12 = C(2,3);
C03 = C(1,4);
C40 = C(5,1);
C31 = C(4,2);
C22 = C(3,3);
C13 = C(2,4); 
C04 = C(1,5);

% compute standardized moments
S00 = 1;
S10 = 0;
S01 = 0;
S20 = 1;
S11 = C11/sqrt(C20*C02);
S02 = 1;
S30 = C30/C20^(3/2);
S21 = C21/C20/sqrt(C02);
S12 = C12/C02/sqrt(C20);
S03 = C03/C02^(3/2);
S40 = C40/C20^2;
S31 = C31/C20^(3/2)/sqrt(C02);
S22 = C22/C20/C02;
S13 = C13/sqrt(C20)/C02^(3/2);
S04 = C04/C02^2;

C4(1) = C00;
C4(2) = C10;
C4(6) = C01;
C4(3) = C20;
C4(7) = C11;
C4(10)= C02;
C4(4) = C30;
C4(8) = C21;
C4(11)= C12;
C4(13)= C03;
C4(5) = C40;
C4(9) = C31;
C4(12)= C22;
C4(14)= C13;
C4(15)= C04;

S4(1) = S00;
S4(2) = S10;
S4(6) = S01;
S4(3) = S20;
S4(7) = S11;
S4(10)= S02;
S4(4) = S30;
S4(8) = S21;
S4(11)= S12;
S4(13)= S03;
S4(5) = S40;
S4(9) = S31;
S4(12)= S22;
S4(14)= S13;
S4(15)= S04;
end