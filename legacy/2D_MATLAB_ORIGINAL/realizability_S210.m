function [S210r,S201r] = realizability_S210(S110,S101,S011,S300,S210,S201,H200,beta)
% realizability_S210 checks and corrects realizablity of S210, S201
%
xr = 1;
X = [S210-S110*S300; S201-S101*S300];
D1 = [1-S101^2, S101*S110-S011; S101*S110-S011, 1-S110^2];
U = sqrtm(D1);
V = U*X;
L = max([0 V'*V]);
dD1 = max([0 det(D1)]);
R = H200*dD1;
if R <= 0 || X'*X < 1000*eps
    xr = 0;
elseif L > R
    xr = sqrt(R/L);
end
Vr = beta*xr*V;
Xr = U\Vr;
S210r = Xr(1) + S110*S300;
S201r = Xr(2) + S101*S300;
end