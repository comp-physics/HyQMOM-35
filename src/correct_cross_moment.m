function C_corrected = correct_cross_moment(S4i, S4j, S3i, S3j, Sij, S22, Cii_Cjj, Cij_sq, scale)
% correct_cross_moment Correct cross moment when complex eigenvalues detected
%
% Inputs:
%   S4i, S4j - 4th order standardized moments
%   S3i, S3j - 3rd order standardized moments  
%   Sij - 2nd order cross standardized moment
%   S22 - current 2nd order diagonal standardized moment
%   Cii_Cjj - product of central moments
%   Cij_sq - squared cross central moment
%   scale - scaling factor (sC_i * sC_j)
%
% Output:
%   C_corrected - corrected central moment

A = S4i + S4j + 2 - 4*S3j^2*S4i - 4*S3i^2*S4j + 8*S3j*Sij*S3i*S22 - 2*S22^2;
B = -2 + 2*S3j^2 - 2*S22 + 2*S3i^2 + 4*S3j*Sij*S3i;
s22min = A/B;
s22min = max(s22min, (2 + 6*S3j*Sij*S3i - S3i^2 - S3j^2)/2);
S22r = min(max(S22, s22min), sqrt(max(eps, Cii_Cjj + 2*Cij_sq))/scale);
C_corrected = S22r * scale;

end
