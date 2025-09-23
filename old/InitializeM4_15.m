function M = InitializeM4_15(M00,umean,vmean,C20,C11,C02)
% InitializeM4 computes fourth-order joint Gaussian moments
%   Input:
%       M00 = number density 
%       umean = M10/M00 (mean velocity u)
%       vmean = M01/M00 (mean velocity v)
%   Output:
%       M = 4th-order moments
%         = [M00,M10,M20,M30,M40,M01,M11,M21,M31,M02,M12,M22,M03,M13,M04]

% standardized moments for Maxwellian
S30 = 0;
S21 = 0;
S12 = 0;
S03 = 0;
S40 = 3;
S31 = 0;
S22 = 1;
S13 = 0;
S04 = 3;

% compute central 4th-order moments (uses rotation to initilize C and M)
C4 = S4toC4(C20,C11,C02,S30,S21,S12,S03,S40,S31,S22,S13,S04);

C20 = C4(3,1);
C11 = C4(2,2);
C02 = C4(1,3);
C30 = C4(4,1);
C21 = C4(3,2);
C12 = C4(2,3);
C03 = C4(1,4);
C40 = C4(5,1);
C31 = C4(4,2);
C22 = C4(3,3);
C13 = C4(2,4);
C04 = C4(1,5);

% compute 4th-order moments
M4 = C4toM4(M00,umean,vmean,C20,C11,C02,C30,C21,C12,C03,C40,C31,C22,C13,C04);

M00 = M4(1,1);
M10 = M4(2,1);
M01 = M4(1,2);
M20 = M4(3,1);
M11 = M4(2,2);
M02 = M4(1,3);
M30 = M4(4,1);
M21 = M4(3,2);
M12 = M4(2,3);
M03 = M4(1,4);
M40 = M4(5,1);
M31 = M4(4,2);
M22 = M4(3,3);
M13 = M4(2,4);
M04 = M4(1,5);

% M = [M00,M10,M20,M30,M40,M01,M11,M21,M31,M02,M12,M22,M03,M13,M04]
M = zeros(15,1);
M(1) = M00;
M(2) = M10;
M(3) = M20;
M(4) = M30;
M(5) = M40;
M(6) = M01;
M(7) = M11;
M(8) = M21;
M(9) = M31;
M(10)= M02;
M(11)= M12;
M(12)= M22;
M(13)= M03;
M(14)= M13;
M(15)= M04;
end