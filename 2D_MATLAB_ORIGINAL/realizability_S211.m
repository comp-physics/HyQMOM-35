function S211r = realizability_S211(e11,e22,e33,e12,e13,d23,S211,beta)
% realizability_S211
%
S211r = S211;
b211 = e12*e13;
G211 = max([0 (e11*e22-e13^2)*(e11*e33-e12^2)]);
sG211 = beta*sqrt(G211);
s211min = d23+(b211-sG211)/e11;
s211max = d23+(b211+sG211)/e11;
%
if S211 <= s211min
    S211r = s211min;
elseif S211 >= s211max
    S211r = s211max;
end
%
end