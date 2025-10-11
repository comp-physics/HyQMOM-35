function [S220r] = realizability_S310_220(S110,S101,S011,S210,S120,S111,S220)
% realizability_S310_220 checks and corrects realizablity of S220 for S310
%
S220r = S220;
D1 = [1-S101^2, S101*S110-S011; S101*S110-S011, 1-S110^2];
dD1 = det(D1);
%
D2 = [1-S011^2, S011*S101-S110, S011*S110-S101;...
      S011*S101-S110, 1-S101^2, S101*S110-S011;...
      S011*S110-S101, S101*S110-S011, 1-S110^2];
V2 = [S210; S120; S111];
L2 = V2'*D2*V2/dD1;
S220 = S110^2 + L2 + 1000*eps;
%
if S220 > S220r
    S220r = S220;
end
%
end