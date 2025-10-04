function [S110, S101, S011, S300, S030, S003, S400, S040, S004, ...
          S210, S201, S120, S111, S102, S021, S012, ...
          S310, S301, S220, S211, S202, S130, S121, S112, S103, S031, S022, S013] = ...
    edge_corner_correction(R110, R101, R011, ...
                           S110r, S101r, S011r, ...
                           S300r1, S030r1, S400r1, S040r1, ...
                           S300r2, S003r2, S400r2, S004r2, ...
                           S030r3, S003r3, S040r3, S004r3, ...
                           S300r, S030r, S003r, S400r, S040r, S004r, ...
                           S210r, S201r, S120r, S111r, S102r, S021r, S012r, ...
                           S310r, S301r, S220r, S211r, S202r, S130r, S121r, S112r, S103r, ...
                           S031r, S022r, S013r)
%EDGE_CORNER_CORRECTION Corrects moments at edges/corners of realizability domain
%
%   Handles cases where one or more 2D correlations (R110, R101, R011) are
%   non-realizable (≤ 0), placing the state at edges or corners of the
%   realizability domain.
%
%   Refactored from 498-line script with ~70% duplication to 243-line function.
%
%   Inputs:
%       R110, R101, R011 - Realizability indicators (1 - S###^2)
%       S###r            - Corrected moments from check2D_all_planes
%
%   Outputs:
%       S###             - Edge/corner-corrected standardized moments

% Route to appropriate edge or corner correction
if R110 <= 0 && R101 > 0 && R011 > 0
    % Edge 1: S110 = ±1 (xy-plane boundary)
    S110 = sign(S110r);
    Smean = (S011r + S101r)/2;
    S011 = sign(S011r) * Smean;
    S101 = sign(S101r) * Smean;
    [S110, S101, S011, ~] = realizability('S2', S110, S101, S011);
    
    S300 = S300r1; S030 = S030r1; S400 = S400r1; S040 = S040r1;
    S003 = S003r; S004 = S004r;
    
    [S210, S201, S120, S111, S102, S021, S012, S310, S301, S220, S211, S202, S130, S121, S112, S103, S031, S022, S013] = ...
        average_edge_variants_xy(S110, S030, S040, S300r, S400r, S021r, S012r, S013r, S031r, S022r, ...
                                  S201r, S102r, S103r, S301r, S202r, S003r, S004r);
    
elseif R101 <= 0 && R110 > 0 && R011 > 0
    % Edge 2: S101 = ±1 (xz-plane boundary)
    S101 = sign(S101r);
    Smean = (S011r + S110r)/2;
    S011 = sign(S011r) * Smean;
    S110 = sign(S110r) * Smean;
    [S110, S101, S011, ~] = realizability('S2', S110, S101, S011);
    
    S300 = S300r2; S003 = S003r2; S400 = S400r2; S004 = S004r2;
    S030 = S030r; S040 = S040r;
    
    [S210, S201, S120, S111, S102, S021, S012, S310, S301, S220, S211, S202, S130, S121, S112, S103, S031, S022, S013] = ...
        average_edge_variants_xz(S101, S003, S004, S300r, S400r, S021r, S012r, S013r, S031r, S022r, ...
                                  S030r, S210r, S120r, S220r, S130r, S310r);
    
elseif R011 <= 0 && R101 > 0 && R110 > 0
    % Edge 3: S011 = ±1 (yz-plane boundary)
    S011 = sign(S011r);
    Smean = (S101r + S110r)/2;
    S101 = sign(S101r) * Smean;
    S110 = sign(S110r) * Smean;
    [S110, S101, S011, ~] = realizability('S2', S110, S101, S011);
    
    S030 = S030r3; S003 = S003r3; S040 = S040r3; S004 = S004r3;
    S300 = S300r; S400 = S400r;
    
    [S210, S201, S120, S111, S102, S021, S012, S310, S301, S220, S211, S202, S130, S121, S112, S103, S031, S022, S013] = ...
        average_edge_variants_yz(S011, S003, S004, S030r, S040r, S201r, S102r, S103r, S301r, S202r, ...
                                  S300r, S210r, S120r, S220r, S130r, S310r);
    
else
    % Corner: all three correlations at boundaries
    S110 = sign(S110r);
    S101 = sign(S101r);
    S011 = sign(S011r);
    if S011 * S101 * S110 ~= 1
        warning('edge_corner_correction: S011*S101*S110 ~= 1 at corner');
    end
    [S110, S101, S011, ~] = realizability('S2', S110, S101, S011);
    
    % Average univariates from three planes
    S300 = (S030r1 + S030r + S300r)/3;
    S030 = (S030r1 + S030r + S030r3)/3;
    S003 = (S003r + S003r2 + S003r3)/3;
    S400 = (S040r1 + S040r + S400r)/3;
    S040 = (S040r1 + S040r + S040r3)/3;
    S004 = (S004r + S004r2 + S004r3)/3;
    
    [S210, S201, S120, S111, S102, S021, S012, S310, S301, S220, S211, S202, S130, S121, S112, S103, S031, S022, S013] = ...
        average_corner_variants(S110, S101, S011, S030r1, S040r1, S030r, S040r, S300r, S400r, S003r3, S004r3);
end

% Validate realizability
S2 = 1 + 2*S110*S101*S011 - (S110^2 + S101^2 + S011^2);
if S2 < 0
    warning('edge_corner_correction: S2 < 0 = %g', S2);
end

end

%% Helper: Edge variant averaging for xy-plane (S110 edge)
function [S210, S201, S120, S111, S102, S021, S012, S310, S301, S220, S211, S202, S130, S121, S112, S103, S031, S022, S013] = ...
    average_edge_variants_xy(S110, S030, S040, S300, S400, S021r, S012r, S013r, S031r, S022r, ...
                              S201r, S102r, S103r, S301r, S202r, S003r, S004r)
% Variant a: SIJK = S110^I * S0(J+I)K
% Variant b: SIJK = S110^J * S(I+J)0K

S210 = (S030 + S110*S300)/2;
S201 = (S021r + S201r)/2;
S120 = (S110*S030 + S300)/2;
S111 = (S110*S021r + S110*S201r)/2;
S102 = (S110*S012r + S102r)/2;
S021 = (S021r + S201r)/2;
S012 = (S012r + S110*S102r)/2;

S310 = (S110*S040 + S110*S400)/2;
S301 = (S110*S031r + S301r)/2;
S220 = (S040 + S400)/2;
S211 = (S031r + S110*S301r)/2;
S202 = (S022r + S202r)/2;
S130 = (S110*S040 + S110*S400)/2;
S121 = (S110*S031r + S301r)/2;
S112 = (S110*S022r + S110*S202r)/2;
S103 = (S110*S013r + S103r)/2;
S031 = (S031r + S110*S301r)/2;
S022 = (S022r + S202r)/2;
S013 = (S013r + S110*S103r)/2;
end

%% Helper: Edge variant averaging for xz-plane (S101 edge)
function [S210, S201, S120, S111, S102, S021, S012, S310, S301, S220, S211, S202, S130, S121, S112, S103, S031, S022, S013] = ...
    average_edge_variants_xz(S101, S003, S004, S300, S400, S021r, S012r, S013r, S031r, S022r, ...
                              S030r, S210r, S120r, S220r, S130r, S310r)
% Variant a: SIJK = S101^I * S0J(K+I)
% Variant b: SIJK = S101^K * S(I+K)J0

S210 = (S013r + S210r)/2;
S201 = (S003 + S101*S300)/2;
S120 = (S101*S021r + S120r)/2;
S111 = (S101*S012r + S101*S210r)/2;
S102 = (S101*S003 + S300)/2;
S021 = (S021r + S101*S120r)/2;
S012 = (S012r + S210r)/2;

S310 = (S101*S013r + S310r)/2;
S301 = (S101*S004 + S101*S400)/2;
S220 = (S022r + S220r)/2;
S211 = (S013r + S101*S310r)/2;
S202 = (S004 + S400)/2;
S130 = (S101*S031r + S130r)/2;
S121 = (S101*S022r + S101*S220r)/2;
S112 = (S101*S013r + S310r)/2;
S103 = (S101*S004 + S101*S400)/2;
S031 = (S031r + S101*S130r)/2;
S022 = (S022r + S220r)/2;
S013 = (S013r + S101*S310r)/2;
end

%% Helper: Edge variant averaging for yz-plane (S011 edge)
function [S210, S201, S120, S111, S102, S021, S012, S310, S301, S220, S211, S202, S130, S121, S112, S103, S031, S022, S013] = ...
    average_edge_variants_yz(S011, S003, S004, S030, S040, S201r, S102r, S103r, S301r, S202r, ...
                              S300r, S210r, S120r, S220r, S130r, S310r)
% Variant a: SIJK = S011^J * SI0(K+J)
% Variant b: SIJK = S011^K * SI(J+K)0

S210 = (S011*S201r + S210r)/2;
S201 = (S201r + S011*S210r)/2;
S120 = (S102r + S120r)/2;
S111 = (S011*S102r + S011*S120r)/2;
S102 = (S102r + S120r)/2;
S021 = (S003 + S011*S030)/2;
S012 = (S011*S003 + S030)/2;

S310 = (S011*S301r + S310r)/2;
S301 = (S301r + S011*S310r)/2;
S220 = (S202r + S220r)/2;
S211 = (S011*S202r + S011*S220r)/2;
S202 = (S202r + S220r)/2;
S130 = (S011*S103r + S130r)/2;
S121 = (S103r + S011*S130r)/2;
S112 = (S011*S103r + S130r)/2;
S103 = (S103r + S011*S130r)/2;
S031 = (S011*S004 + S011*S040)/2;
S022 = (S004 + S040)/2;
S013 = (S011*S004 + S011*S040)/2;
end

%% Helper: Corner variant averaging (all three boundaries)
function [S210, S201, S120, S111, S102, S021, S012, S310, S301, S220, S211, S202, S130, S121, S112, S103, S031, S022, S013] = ...
    average_corner_variants(S110, S101, S011, S030a, S040a, S030b, S040b, S300c, S400c, S003c, S004c)
% Variant a: SIJK = S110^J * S101^K * S(I+J+K)00
% Variant b: SIJK = S110^I * S011^K * S0(I+J+K)0
% Variant c: SIJK = S101^I * S011^J * S00(I+J+K)

% Variant a
S210a = S110*S030a; S201a = S101*S030a; S120a = S030a; S111a = S110*S101*S030a; S102a = S030a;
S021a = S101*S030a; S012a = S110*S030a;
S310a = S110*S040a; S301a = S101*S040a; S220a = S040a; S211a = S110*S101*S040a; S202a = S040a;
S130a = S110*S040a; S121a = S101*S040a; S112a = S110*S040a; S103a = S101*S040a;
S031a = S110*S101*S040a; S022a = S040a; S013a = S110*S101*S040a;

% Variant b
S210b = S030b; S201b = S011*S030b; S120b = S110*S030b; S111b = S110*S011*S030b; S102b = S110*S030b;
S021b = S011*S030b; S012b = S030b;
S310b = S110*S040b; S301b = S110*S011*S040b; S220b = S040b; S211b = S011*S040b; S202b = S040b;
S130b = S110*S040b; S121b = S110*S011*S040b; S112b = S110*S040b; S103b = S110*S011*S040b;
S031b = S011*S040b; S022b = S040b; S013b = S011*S040b;

% Variant c
S210c = S011*S003c; S201c = S003c; S120c = S101*S003c; S111c = S101*S011*S003c; S102c = S101*S003c;
S021c = S003c; S012c = S011*S003c;
S310c = S101*S011*S004c; S301c = S101*S004c; S220c = S004c; S211c = S011*S004c; S202c = S004c;
S130c = S101*S011*S004c; S121c = S101*S004c; S112c = S101*S011*S004c; S103c = S101*S004c;
S031c = S011*S004c; S022c = S004c; S013c = S011*S004c;

% Average all three
S210 = (S210a + S210b + S210c)/3;
S201 = (S201a + S201b + S201c)/3;
S120 = (S120a + S120b + S120c)/3;
S111 = (S111a + S111b + S111c)/3;
S102 = (S102a + S102b + S102c)/3;
S021 = (S021a + S021b + S021c)/3;
S012 = (S012a + S012b + S012c)/3;

S310 = (S310a + S310b + S310c)/3;
S301 = (S301a + S301b + S301c)/3;
S220 = (S220a + S220b + S220c)/3;
S211 = (S211a + S211b + S211c)/3;
S202 = (S202a + S202b + S202c)/3;
S130 = (S130a + S130b + S130c)/3;
S121 = (S121a + S121b + S121c)/3;
S112 = (S112a + S112b + S112c)/3;
S103 = (S103a + S103b + S103c)/3;

S031 = (S031a + S031b + S031c)/3;
S022 = (S022a + S022b + S022c)/3;
S013 = (S013a + S013b + S013c)/3;
end
