function S220r = realizablity_S220(S110,S220,A220)
% realizablity_S220 checks maximum bounds and corrects S220
%
S220r = S220;
s220min = max([S110^2 1-A220]);
s220max = 1+A220;
if S220 < s220min
    S220r = s220min;
elseif S220 > s220max
    S220r = s220max;
end
end