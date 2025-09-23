function [M] = Moment5_nr(M4)
% Moment5 computes 2-D Grad-HyQMOM closure for 5th-order moments w/o
% rotation
%   Input:
%       M4 = all moments up to 4th order 
%          = [M00,M10,M20,M30,M40,M01,M11,M21,M31,M02,M12,M22,M03,M13,M04]
%   Output:
%       M = 5th-order moments found from closure
%         = [M50 M41 M32 M23 M14 M05]

% M00 = M4(1);
% M10 = M4(2);
% M01 = M4(3);
% M20 = M4(4);
% M11 = M4(5);
% M02 = M4(6);
% M30 = M4(7);
% M21 = M4(8);
% M12 = M4(9);
% M03 = M4(10);
% M40 = M4(11);
% M31 = M4(12);
% M22 = M4(13);
% M13 = M4(14);
% M04 = M4(15);

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
C4 = M4toC4(M00,M10,M01,M20,M11,M02,M30,M21,M12,M03,M40,M31,M22,M13,M04);

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

% compute standardized moments
S4 = C4toS4_nr(C20,C11,C02,C30,C21,C12,C03,C40,C31,C22,C13,C04);

S30 = S4(4,1);
S21 = S4(3,2);
S12 = S4(2,3);
S03 = S4(1,4);
S40 = S4(5,1);
S31 = S4(4,2);
S22 = S4(3,3);
S13 = S4(2,4);
S04 = S4(1,5);

% compute standardized 5th-order moments using 2-D Grad-HyQMOM
S5 = S5flux_nr(S30,S21,S12,S03,S40,S31,S22,S13,S04);

S50 = S5(1);
S41 = S5(2); 
S32 = S5(3); 
S23 = S5(4); 
S14 = S5(5);
S05 = S5(6);

% compute central 5th-order moments
C5 = S5toC5_nr(C20,C11,C02,S30,S21,S12,S03,S40,S31,S22,S13,S04,S50,S41,S32,S23,S14,S05);

C50 = C5(6,1);
C41 = C5(5,2); 
C32 = C5(4,3); 
C23 = C5(3,4); 
C14 = C5(2,5);
C05 = C5(1,6);

% compute 5th-order moments
M5 = C5toM5(M00,M10/M00,M01/M00,C20,C11,C02,C30,C21,C12,C03,C40,C31,C22,C13,C04,C50,C41,C32,C23,C14,C05);

M50 = M5(6,1);
M41 = M5(5,2); 
M32 = M5(4,3); 
M23 = M5(3,4); 
M14 = M5(2,5);
M05 = M5(1,6);

M = [M50 M41 M32 M23 M14 M05];

end