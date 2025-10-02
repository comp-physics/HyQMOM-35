function [S310r,S220r] = realizability_S310(S110,S101,S011,S300,S210,S201,S120,S111,S310,S220,H200,beta)
% realizability_S310 checks and corrects realizablity of S310 and S220
%
S310r = S310;
S220r = S220;
%
b310 = S111*((1 - S110^2)*S201 + (S101*S110 - S011)*S210 + (S011*S110 - S101)*S300) ...
     + S120*((1 - S101^2)*S210 + (S011*S101 - S110)*S300 + (S101*S110 - S011)*S201) ...
     + S210*((1 - S011^2)*S300 + (S011*S101 - S110)*S210 + (S011*S110 - S101)*S201);
%
D1 = [1-S101^2, S101*S110-S011;S101*S110-S011, 1-S110^2];
dD1 = det(D1);
V1 = [S210-S110*S300; S201-S101*S300];
L1 = V1'*D1*V1/dD1;
%
G310b = H200 - L1;
if G310b < 0
    warning('G310b < 0')
    G310b = 0;
end
%
D2 = [1-S011^2, S011*S101-S110, S011*S110-S101;...
      S011*S101-S110, 1-S101^2, S101*S110-S011;...
      S011*S110-S101, S101*S110-S011, 1-S110^2];
V2 = [S210; S120; S111];
L2 = V2'*D2*V2/dD1;
%
G310a = S220 - S110^2 - L2 + 1000*eps;
if G310a < 0 || L2 < 0
    warning('G310a < 0 || L2 < 0')
    % [G310a L2 b310/dD1]
    % S220r = S110^2 + L2 + 1000*eps;
    G310a = 0;
end
%
G310 = G310a*G310b;
sG310 = beta*sqrt(G310);
%
s310min = S110 + b310/dD1 - sG310;
s310max = S110 + b310/dD1 + sG310;
%
if S310 <= s310min
    S310r = s310min;
elseif S310 >= s310max
    S310r = s310max;
end
%
end