function [S30,S40,S11,S21,S31,S12,S22,S03,S13,S04] = check2D(S30,S40,S11,S21,S31,S12,S22,S03,S13,S04)
% check2D checks and corrects 2D moments in 3D code
%
h2min = 1000*eps;
if abs(S11) >= 1
    % Collapse both S11 >= 1 and S11 <= -1 branches
    S11 = sign(S11);
    s3m = sqrt(abs(S30*S03));
    s4m = sqrt(S40*S04);
    H2m = max(h2min, s4m - s3m^2 - 1);
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
    % check and correct realizability of cross moments
    [S21,S12,S31,S22,S13] = realizable_2D(S30,S40,S11,S21,S31,S12,S22,S03,S13,S04);
end
end