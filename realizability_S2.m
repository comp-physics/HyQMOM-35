function [S110r,S101r,S011r,S2r] = realizability_S2(S110,S101,S011)
% realizability_S2 checks and corrects realizablity of 2nd-order moments
%
S2 = 1 + 2*S110*S101*S011 - (S110^2+S101^2+S011^2);
xr = 1;
if S2 < 0
    Y = @(x) 1 + 2*S110*S101*S011*x.^3 - (S110^2+S101^2+S011^2)*x.^2;
    xr = fzero(Y,[0 1]);
end
xr = 0.9999*xr;
S110r = xr*S110;
S101r = xr*S101;
S011r = xr*S011;
S2r = 1 + 2*S110r*S101r*S011r - (S110r^2+S101r^2+S011r^2);
if S2r < 0
    warning('S2 < 0 after correction in realizability_S2')
    disp([S2,S2r])
end
end