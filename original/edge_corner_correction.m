% treat cases where one or more 2D 2nd-order moments is non-realizable (corners and edges)
%
if R110 <= 0 && R101 > 0 && R011 > 0 % treat edge with S110 = \pm 1
    S110r = sign(S110r);
    % R011 = R101 on edge
    Smean = (S011r+S101r)/2;
    S011r = sign(S011r)*Smean;
    S101r = sign(S101r)*Smean;
    % check and correct realizability of S110, S101, S011
    [S110r,S101r,S011r,S2] = realizability_S2(S110r,S101r,S011r);
    %
    S300r = S300r1;
    S030r = S030r1;
    S400r = S400r1;
    S040r = S040r1;
    %
    % SIJK = S110^I*S0(J+I)K;
    % S110a = S110r;
    % S101a = S110r*S011r;
    % S011a =       S011r;
    %
    S300a = S110r*S030r;
    S210a =       S030r;
    S201a =       S021r;
    S120a = S110r*S030r;
    S111a = S110r*S021r;
    S102a = S110r*S012r;
    S030a =       S030r;
    S021a =       S021r;
    S012a =       S012r;
    S003a =       S003r;
    %
    S400a =       S040r;
    S310a = S110r*S040r;
    S301a = S110r*S031r;
    S220a =       S040r;
    S211a =       S031r;
    S202a =       S022r;
    S130a = S110r*S040r;
    S121a = S110r*S031r;
    S112a = S110r*S022r;
    S103a = S110r*S013r;
    S040a =       S040r;
    S031a =       S031r;
    S022a =       S022r;
    S013a =       S013r;
    S004a =       S004r;
    %
    % SIJK = S110^J*S(I+J)0K;
    S110b = S110r;
    S101b =       S101r;
    S011b = S110r*S101r;
    %
    S300b =       S300r;
    S210b = S110r*S300r;
    S201b =       S201r;
    S120b =       S300r;
    S111b = S110r*S201r;
    S102b =       S102r;
    S030b = S110r*S300r;
    S021b =       S201r;
    S012b = S110r*S102r;
    S003b =       S003r;
    %
    S400b =       S400r;
    S310b = S110r*S400r;
    S301b =       S301r;
    S220b =       S400r;
    S211b = S110r*S301r;
    S202b =       S202r;
    S130b = S110r*S400r;
    S121b =       S301r;
    S112b = S110r*S202r;
    S103b =       S103r;
    S040b =       S400r;
    S031b = S110r*S301r;
    S022b =       S202r;
    S013b = S110r*S103r;
    S004b =       S004r;
    %
    % S110r = (S110a+S110b)/2;
    % S101r = (S101a+S101b)/2;
    % S011r = (S011a+S011b)/2;
    %
    S300r = (S300a+S300b)/2;
    S210r = (S210a+S210b)/2;
    S201r = (S201a+S201b)/2;
    S120r = (S120a+S120b)/2;
    S111r = (S111a+S111b)/2;
    S102r = (S102a+S102b)/2;
    S030r = (S030a+S030b)/2;
    S021r = (S021a+S021b)/2;
    S012r = (S012a+S012b)/2;
    S003r = (S003a+S003b)/2;
    %
    S400r = (S400a+S400b)/2;
    S310r = (S310a+S310b)/2;
    S301r = (S301a+S301b)/2;
    S220r = (S220a+S220b)/2;
    S211r = (S211a+S211b)/2;
    S202r = (S202a+S202b)/2;
    S130r = (S130a+S130b)/2;
    S121r = (S121a+S121b)/2;
    S112r = (S112a+S112b)/2;
    S103r = (S103a+S103b)/2;
    S040r = (S040a+S040b)/2;
    S031r = (S031a+S031b)/2;
    S022r = (S022a+S022b)/2;
    S013r = (S013a+S013b)/2;
    S004r = (S004a+S004b)/2;
    %
elseif R101 <= 0 && R110 > 0 && R011 > 0 % treat edge with S101 = \pm 1
    S101r = sign(S101r);
    % R011 = R110 on edge
    Smean = (S011r+S110r)/2;
    S011r = sign(S011r)*Smean;
    S110r = sign(S110r)*Smean;
    % check and correct realizability of S110, S101, S011
    [S110r,S101r,S011r,S2] = realizability_S2(S110r,S101r,S011r);
    %
    S300r = S300r2;
    S003r = S003r2;
    S400r = S400r2;
    S004r = S004r2;
    %
    % SIJK = S101^I*S0J(K+I);
    % S110a = S101r*S011r;
    % S101a = S101r;
    % S011a =       S011r;
    %
    S300a = S101r*S003r;
    S210a =       S013r;
    S201a =       S003r;
    S120a = S101r*S021r;
    S111a = S101r*S012r;
    S102a = S101r*S003r;
    S030a =       S030r;
    S021a =       S021r;
    S012a =       S012r;
    S003a =       S003r;
    %
    S400a =       S004r;
    S310a = S101r*S013r;
    S301a = S101r*S004r;
    S220a =       S022r;
    S211a =       S013r;
    S202a =       S004r;
    S130a = S101r*S031r;
    S121a = S101r*S022r;
    S112a = S101r*S013r;
    S103a = S101r*S004r;
    S040a =       S040r;
    S031a =       S031r;
    S022a =       S022r;
    S013a =       S013r;
    S004a =       S004r;
    %
    % SIJK = S101^K*S(I+K)J0;
    % S110b =       S110r;
    % S101b = S101r;
    % S011b = S101r*S110r;
    %
    S300b =       S300r;
    S210b =       S210r;
    S201b = S101r*S300r;
    S120b =       S120r;
    S111b = S101r*S210r;
    S102b =       S300r;
    S030b =       S030r;
    S021b = S101r*S120r;
    S012b =       S210r;
    S003b = S101r*S300r;
    %
    S400b =       S400r;
    S310b =       S310r;
    S301b = S101r*S400r;
    S220b =       S220r;
    S211b = S101r*S310r;
    S202b =       S400r;
    S130b =       S130r;
    S121b = S101r*S220r;
    S112b =       S310r;
    S103b = S101r*S400r;
    S040b =       S040r;
    S031b = S101r*S130r;
    S022b =       S220r;
    S013b = S101r*S310r;
    S004b =       S400r;
    %
    % S110r = (S110a+S110b)/2;
    % S101r = (S101a+S101b)/2;
    % S011r = (S011a+S011b)/2;
    %
    S300r = (S300a+S300b)/2;
    S210r = (S210a+S210b)/2;
    S201r = (S201a+S201b)/2;
    S120r = (S120a+S120b)/2;
    S111r = (S111a+S111b)/2;
    S102r = (S102a+S102b)/2;
    S030r = (S030a+S030b)/2;
    S021r = (S021a+S021b)/2;
    S012r = (S012a+S012b)/2;
    S003r = (S003a+S003b)/2;
    %
    S400r = (S400a+S400b)/2;
    S310r = (S310a+S310b)/2;
    S301r = (S301a+S301b)/2;
    S220r = (S220a+S220b)/2;
    S211r = (S211a+S211b)/2;
    S202r = (S202a+S202b)/2;
    S130r = (S130a+S130b)/2;
    S121r = (S121a+S121b)/2;
    S112r = (S112a+S112b)/2;
    S103r = (S103a+S103b)/2;
    S040r = (S040a+S040b)/2;
    S031r = (S031a+S031b)/2;
    S022r = (S022a+S022b)/2;
    S013r = (S013a+S013b)/2;
    S004r = (S004a+S004b)/2;
    %
elseif R011 <= 0 && R101 > 0 && R110 > 0 % treat edge with S011 = \pm 1
    S011r = sign(S011r);
    % R101 = R110 on edge
    Smean = (S101r+S110r)/2;
    S101r = sign(S101r)*Smean;
    S110r = sign(S110r)*Smean;
    % check and correct realizability of S110, S101, S011
    [S110r,S101r,S011r,S2] = realizability_S2(S110r,S101r,S011r);
    %
    S030r = S030r3;
    S003r = S003r3;
    S040r = S040r3;
    S004r = S004r3;
    %
    % SIJK = S011^J*SI0(K+J);
    % S110a = S011r*S101r;
    % S101a =       S101r;
    % S011a = S011r;
    %
    S300a =       S300r;
    S210a = S011r*S201r;
    S201a =       S201r;
    S120a =       S102r;
    S111a = S011r*S102r;
    S102a =       S102r;
    S030a = S011r*S003r;
    S021a =       S003r;
    S012a = S011r*S003r;
    S003a =       S003r;
    %
    S400a =       S400r;
    S310a = S011r*S301r;
    S301a =       S301r;
    S220a =       S202r;
    S211a = S011r*S202r;
    S202a =       S202r;
    S130a = S011r*S103r;
    S121a =       S103r;
    S112a = S011r*S103r;
    S103a =       S103r;
    S040a =       S004r;
    S031a = S011r*S004r;
    S022a =       S004r;
    S013a = S011r*S004r;
    S004a =       S004r;
    %
    % SIJK = S011^K*SI(J+K)0;
    % S110b =       S110r;
    % S101b = S011r*S110r;
    % S011b = S011r;
    %
    S300b =       S300r;
    S210b =       S210r;
    S201b = S011r*S210r;
    S120b =       S120r;
    S111b = S011r*S120r;
    S102b =       S120r;
    S030b =       S030r;
    S021b = S011r*S030r;
    S012b =       S030r;
    S003b = S011r*S030r;
    %
    S400b =       S400r;
    S310b =       S310r;
    S301b = S011r*S310r;
    S220b =       S220r;
    S211b = S011r*S220r;
    S202b =       S220r;
    S130b =       S130r;
    S121b = S011r*S130r;
    S112b =       S130r;
    S103b = S011r*S130r;
    S040b =       S040r;
    S031b = S011r*S040r;
    S022b =       S040r;
    S013b = S011r*S040r;
    S004b =       S040r;
    %
    % S110r = (S110a+S110b)/2;
    % S101r = (S101a+S101b)/2;
    % S011r = (S011a+S011b)/2;
    %
    S300r = (S300a+S300b)/2;
    S210r = (S210a+S210b)/2;
    S201r = (S201a+S201b)/2;
    S120r = (S120a+S120b)/2;
    S111r = (S111a+S111b)/2;
    S102r = (S102a+S102b)/2;
    S030r = (S030a+S030b)/2;
    S021r = (S021a+S021b)/2;
    S012r = (S012a+S012b)/2;
    S003r = (S003a+S003b)/2;
    %
    S400r = (S400a+S400b)/2;
    S310r = (S310a+S310b)/2;
    S301r = (S301a+S301b)/2;
    S220r = (S220a+S220b)/2;
    S211r = (S211a+S211b)/2;
    S202r = (S202a+S202b)/2;
    S130r = (S130a+S130b)/2;
    S121r = (S121a+S121b)/2;
    S112r = (S112a+S112b)/2;
    S103r = (S103a+S103b)/2;
    S040r = (S040a+S040b)/2;
    S031r = (S031a+S031b)/2;
    S022r = (S022a+S022b)/2;
    S013r = (S013a+S013b)/2;
    S004r = (S004a+S004b)/2;
    %
else % treat corner
    S110r = sign(S110r);
    S101r = sign(S101r);
    S011r = sign(S011r);
    if S011r*S101r*S110r ~= 1
        warning('S011r*S101r*S110r ~= 1 at corner')
    end
    % check and correct realizability of S110, S101, S011
    [S110r,S101r,S011r,S2] = realizability_S2(S110r,S101r,S011r);
    %
    S300a = S030r1;
    S300b = S030r2;
    S300c = S300r;
    S030a = S030r1;
    S030b = S030r;
    S030c = S030r3;
    S003a = S003r;
    S003b = S003r2;
    S003c = S003r3;
    %
    S400a = S040r1;
    S400b = S040r2;
    S400c = S040r;
    S040a = S040r1;
    S040b = S040r;
    S040c = S040r3;
    S004a = S004r;
    S004b = S004r2;
    S004c = S004r3;
    %
    % SIJK = S110^J*S101^K*S(I+J+K)00;
    S210a = S110r      *S300a;
    S120a =             S300a;
    S201a =       S101r*S300a;
    S102a =             S300a;
    S021a =       S101r*S300a;
    S012a = S110r      *S300a;
    %
    S111a = S110r*S101r*S300a;
    %
    S310a = S110r      *S400a;
    S130a = S110r      *S400a;
    S103a =       S101r*S400a;
    S301a =       S101r*S400a;
    S031a = S110r*S101r*S400a;
    S013a = S110r*S101r*S400a;
    %
    S220a =             S400a;
    S022a =             S400a;
    S202a =             S400a;
    %
    S211a = S110r*S101r*S400a;
    S121a =       S101r*S400a;
    S112a = S110r      *S400a;
    %
    % SIJK = S110^I*S011^K*S0(I+J+K)0;
    S210b =             S030b;
    S120b = S110r      *S030b;
    S201b =       S011r*S030b;
    S102b = S110r      *S030b;
    S021b =       S011r*S030b;
    S012b =             S030b;
    %
    S111b = S110r*S011r*S030b;
    %
    S310b = S110r      *S040b;
    S130b = S110r      *S040b;
    S103b = S110r*S011r*S040b;
    S301b = S110r*S011r*S040b;
    S031b =       S011r*S040b;
    S013b =       S011r*S040b;
    %
    S220b =             S040b;
    S022b =             S040b;
    S202b =             S040b;
    %
    S211b =       S011r*S040b;
    S121b = S110r*S011r*S040b;
    S112b = S110r      *S040b;
    %
    % SIJK = S101^I*S011^J*S00(I+J+K);
    S210c =       S011r*S003c;
    S120c = S101r      *S003c;
    S201c =             S003c;
    S102c = S101r      *S003c;
    S021c =             S003c;
    S012c =       S011r*S003c;
    %
    S111c = S101r*S011r*S003c;
    %
    S310c = S101r*S011r*S004c;
    S130c = S101r*S011r*S004c;
    S103c = S101r      *S004c;
    S301c = S101r      *S004c;
    S031c =       S011r*S004c;
    S013c =       S011r*S004c;
    %
    S220c =             S004c;
    S022c =             S004c;
    S202c =             S004c;
    %
    S211c =       S011r*S004c;
    S121c = S101r      *S004c;
    S112c = S101r*S011r*S004c;
    %
    % average
    S300r = (S300a+S300b+S300c)/3;
    S210r = (S210a+S210b+S210c)/3;
    S201r = (S201a+S201b+S201c)/3;
    S120r = (S120a+S120b+S120c)/3;
    S111r = (S111a+S111b+S111c)/3;
    S102r = (S102a+S102b+S102c)/3;
    S030r = (S030a+S030b+S030c)/3;
    S021r = (S021a+S021b+S021c)/3;
    S012r = (S012a+S012b+S012c)/3;
    S003r = (S003a+S003b+S003c)/3;
    %
    S400r = (S400a+S400b+S400c)/3;
    S310r = (S310a+S310b+S310c)/3;
    S301r = (S301a+S301b+S301c)/3;
    S220r = (S220a+S220b+S220c)/3;
    S211r = (S211a+S211b+S211c)/3;
    S202r = (S202a+S202b+S202c)/3;
    S130r = (S130a+S130b+S130c)/3;
    S121r = (S121a+S121b+S121c)/3;
    S112r = (S112a+S112b+S112c)/3;
    S103r = (S103a+S103b+S103c)/3;
    S040r = (S040a+S040b+S040c)/3;
    S031r = (S031a+S031b+S031c)/3;
    S022r = (S022a+S022b+S022c)/3;
    S013r = (S013a+S013b+S013c)/3;
    S004r = (S004a+S004b+S004c)/3;
    %
end
% update with edge/corner values
S300=S300r;
S400=S400r;
S110=S110r;
S210=S210r;
S310=S310r;
S120=S120r;
S220=S220r;
S030=S030r;
S130=S130r;
S040=S040r;
S101=S101r;
S201=S201r;
S301=S301r;
S102=S102r;
S202=S202r;
S003=S003r;
S103=S103r;
S004=S004r;
S011=S011r;
S111=S111r;
S211=S211r;
S021=S021r;
S121=S121r;
S031=S031r;
S012=S012r;
S112=S112r;
S013=S013r;
S022=S022r;
%
% S2 = 1 + 2*S110*S101*S011 - (S110^2+S101^2+S011^2);
if S2 < 0
    warning('S2 < 0')
    S2
end