function [v6min,v6max,Mr] = eigenvalues6_hyperbolic_3D(M,direction,flag2D,Ma)
% eigenvalues6_hyperbolic_3D - Unified eigenvalue computation for 3D flux Jacobian
%
% Inputs:
%   M         - 35-component moment vector
%   direction - 'x' or 'y' for direction of flux Jacobian
%   flag2D    - flag for 2D case
%   Ma        - Mach number
%
% Outputs:
%   v6min     - minimum eigenvalue
%   v6max     - maximum eigenvalue
%   Mr        - corrected moment vector (if realizability issues)
%
% M = [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,...
%      M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,...
%      M031,M012,M112,M013,M022]
%
Mr = M;

% Extract all moments
m000 = M(1);
m100 = M(2);
m200 = M(3);
m300 = M(4);
m400 = M(5);
m010 = M(6);
m110 = M(7);
m210 = M(8);
m310 = M(9);
m020 = M(10);
m120 = M(11);
m220 = M(12);
m030 = M(13);
m130 = M(14);
m040 = M(15);
m001 = M(16);
m101 = M(17);
m201 = M(18);
m301 = M(19);
m002 = M(20);
m102 = M(21);
m202 = M(22);
m003 = M(23);
m103 = M(24);
m004 = M(25);
m011 = M(26);
m111 = M(27);
m211 = M(28);
m021 = M(29);
m121 = M(30);
m031 = M(31);
m012 = M(32);
m112 = M(33);
m013 = M(34);
m022 = M(35);

% Compute eigenvalues based on direction
if strcmpi(direction, 'x')
    % X-direction: UV moments
    J6 = jacobian6(m000,m010,m020,m030,m040,m100,m110,m120,m130,m200,m210,m220,m300,m310,m400);
    lam6a = eig(J6);
    lam6ar = sort(real(lam6a));
    v6min = lam6ar(1);
    v6max = lam6ar(6);
    
    % X-direction: UW moments
    J6 = jacobian6(m000,m001,m002,m003,m004,m100,m101,m102,m103,m200,m201,m202,m300,m301,m400);
    lam6b = eig(J6);
    lam6br = sort(real(lam6b));
    v6min = min([v6min lam6br(1)]);
    v6max = max([v6max lam6br(6)]);
    
elseif strcmpi(direction, 'y')
    % Y-direction: VU moments
    J6 = jacobian6(m000,m100,m200,m300,m400,m010,m110,m210,m310,m020,m120,m220,m030,m130,m040);
    lam6a = eig(J6);
    lam6ar = sort(real(lam6a));
    v6min = lam6ar(1);
    v6max = lam6ar(6);
    
    % Y-direction: VW moments
    J6 = jacobian6(m000,m001,m002,m003,m004,m010,m011,m012,m013,m020,m021,m022,m030,m031,m040);
    lam6b = eig(J6);
    lam6br = sort(real(lam6b));
    v6min = min([v6min lam6br(1)]);
    v6max = max([v6max lam6br(6)]);
else
    error('eigenvalues6_hyperbolic_3D: direction must be ''x'' or ''y''');
end

% Check for complex eigenvalues and correct if needed
if max(abs(imag(lam6a))) > 1000*eps || max(abs(imag(lam6b))) > 1000*eps
    % Get mean velocities
    M000 = M(1);
    M100 = M(2);
    M010 = M(6);
    M001 = M(16);
    umean = M100/M000;
    vmean = M010/M000;
    wmean = M001/M000;
    
    % Get central and standardized moments
    [C4,S4] = M2CS4_35(M);
    
    C200=C4(3);
    C300=C4(4);
    C400=C4(5);
    C110=C4(7);
    C210=C4(8);
    C310=C4(9);
    C020=C4(10);
    C120=C4(11);
    C220=C4(12);
    C030=C4(13);
    C130=C4(14);
    C040=C4(15);
    C101=C4(17);
    C201=C4(18);
    C301=C4(19);
    C002=C4(20);
    C102=C4(21);
    C202=C4(22);
    C003=C4(23);
    C103=C4(24);
    C004=C4(25);
    C011=C4(26);
    C111=C4(27);
    C211=C4(28);
    C021=C4(29);
    C121=C4(30);
    C031=C4(31);
    C012=C4(32);
    C112=C4(33);
    C013=C4(34);
    C022=C4(35);
    
    S110 = S4(7);
    S300 = S4(4);
    S030 = S4(13);
    S400 = S4(5);
    S220 = S4(12);
    S040 = S4(15);
    S101 = S4(17);
    S202 = S4(22);
    S011 = S4(26);
    S022 = S4(35);
    S003 = S4(23);
    S004 = S4(25);
    
    % Direction-specific corrections
    if strcmpi(direction, 'x')
        if max(abs(imag(lam6a))) > 1000*eps
            % Correction for UV moments
            A220 = S400 + S040 + 2 - 4*S030^2*S400 - 4*S300^2*S040 + 8*S030*S110*S300*S220 - 2*S220^2;
            B220 = -2 + 2*S030^2 - 2*S220 + 2*S300^2 + 4*S030*S110*S300;
            s22min = A220/B220;
            s22min = max(s22min,(2+6*S030*S110*S300-S300^2-S030^2)/2);
            sC200 = sqrt(max(eps,C200));
            sC020 = sqrt(max(eps,C020));
            S220r = min(max(S220,s22min),sqrt(max(eps,C200*C020+2*C110^2))/sC200/sC020);
            C220r = S220r*sC200*sC020;
        end
        
        if max(abs(imag(lam6b))) > 1000*eps
            % Correction for UW moments
            A202 = S400 + S004 + 2 - 4*S003^2*S400 - 4*S300^2*S004 + 8*S003*S101*S300*S202 - 2*S202^2;
            B202 = -2 + 2*S003^2 - 2*S202 + 2*S300^2 + 4*S003*S101*S300;
            s22min = A202/B202;
            s22min = max(s22min,(2+6*S003*S101*S300-S300^2-S003^2)/2);
            sC200 = sqrt(max(eps,C200));
            sC002 = sqrt(max(eps,C002));
            S202r = min(max(S202,s22min),sqrt(max(eps,C200*C002+2*C101^2))/sC200/sC002);
            C202r = S202r*sC200*sC002;
        end
    else % y-direction
        if max(abs(imag(lam6a))) > 1000*eps
            % Correction for VU moments
            A220 = S400 + S040 + 2 - 4*S030^2*S400 - 4*S300^2*S040 + 8*S030*S110*S300*S220 - 2*S220^2;
            B220 = -2 + 2*S030^2 - 2*S220 + 2*S300^2 + 4*S030*S110*S300;
            s22min = A220/B220;
            s22min = max(s22min,(2+6*S030*S110*S300-S300^2-S030^2)/2);
            sC200 = sqrt(max(eps,C200));
            sC020 = sqrt(max(eps,C020));
            S220r = min(max(S220,s22min),sqrt(max(eps,C200*C020+2*C110^2))/sC200/sC020);
            C220r = S220r*sC200*sC020;
        end
        
        if max(abs(imag(lam6b))) > 1000*eps
            % Correction for VW moments
            A022 = S040 + S004 + 2 - 4*S003^2*S040 - 4*S030^2*S004 + 8*S003*S011*S030*S022 - 2*S022^2;
            B022 = -2 + 2*S003^2 - 2*S022 + 2*S030^2 + 4*S003*S011*S030;
            s22min = A022/B022;
            s22min = max(s22min,(2+6*S003*S011*S030-S030^2-S003^2)/2);
            sC020 = sqrt(max(eps,C020));
            sC002 = sqrt(max(eps,C002));
            S022r = min(max(S022,s22min),sqrt(max(eps,C020*C002+2*C011^2))/sC020/sC002);
            C022r = S022r*sC020*sC002;
        end
    end
    
    % Update corrected moments
    if exist('C220r','var')
        C220 = C220r;
    end
    if exist('C202r','var')
        C202 = C202r;
    end
    if exist('C022r','var')
        C022 = C022r;
    end
    
    % Reconstruct moments
    M4 = C4toM4_3D(M000,umean,vmean,wmean,C200,C110,C101,C020,C011,C002,C300,...
        C210,C201,C120,C111,C102,C030,C021,C012,C003,C400,C310,C301,C220,C211,...
        C202,C130,C121,C112,C103,C040,C031,C022,C013,C004);
    
    % Extract corrected moments
    M000 = M4(1,1,1);
    M100 = M4(2,1,1);
    M010 = M4(1,2,1);
    M001 = M4(1,1,2);
    
    M200 = M4(3,1,1);
    M110 = M4(2,2,1);
    M101 = M4(2,1,2);
    M020 = M4(1,3,1);
    M011 = M4(1,2,2);
    M002 = M4(1,1,3);
    
    M300 = M4(4,1,1);
    M210 = M4(3,2,1);
    M201 = M4(3,1,2);
    M120 = M4(2,3,1);
    M111 = M4(2,2,2);
    M102 = M4(2,1,3);
    M030 = M4(1,4,1);
    M021 = M4(1,3,2);
    M012 = M4(1,2,3);
    M003 = M4(1,1,4);
    
    M400 = M4(5,1,1);
    M310 = M4(4,2,1);
    M301 = M4(4,1,2);
    M220 = M4(3,3,1);
    M211 = M4(3,2,2);
    M202 = M4(3,1,3);
    M130 = M4(2,4,1);
    M121 = M4(2,3,2);
    M112 = M4(2,2,3);
    M103 = M4(2,1,4);
    M040 = M4(1,5,1);
    M031 = M4(1,4,2);
    M022 = M4(1,3,3);
    M013 = M4(1,2,4);
    M004 = M4(1,1,5);
    
    % Pack corrected moments
    Mh = [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,...
          M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,...
          M031,M012,M112,M013,M022];
    
    [~,~,~,Mr] = Flux_closure35_and_realizable_3D(Mh,flag2D,Ma);
    
    % Recompute eigenvalues with corrected moments
    m000 = Mr(1);
    m100 = Mr(2);
    m200 = Mr(3);
    m300 = Mr(4);
    m400 = Mr(5);
    m010 = Mr(6);
    m110 = Mr(7);
    m210 = Mr(8);
    m310 = Mr(9);
    m020 = Mr(10);
    m120 = Mr(11);
    m220 = Mr(12);
    m030 = Mr(13);
    m130 = Mr(14);
    m040 = Mr(15);
    m001 = Mr(16);
    m101 = Mr(17);
    m201 = Mr(18);
    m301 = Mr(19);
    m002 = Mr(20);
    m102 = Mr(21);
    m202 = Mr(22);
    m003 = Mr(23);
    m103 = Mr(24);
    m004 = Mr(25);
    m011 = Mr(26);
    m021 = Mr(29);
    m031 = Mr(31);
    m012 = Mr(32);
    m013 = Mr(34);
    m022 = Mr(35);
    
    % Recompute eigenvalues with corrected moments
    if strcmpi(direction, 'x')
        J6 = jacobian6(m000,m010,m020,m030,m040,m100,m110,m120,m130,m200,m210,m220,m300,m310,m400);
        lam6a = eig(J6);
        lam6ar = sort(real(lam6a));
        v6min = lam6ar(1);
        v6max = lam6ar(6);
        
        J6 = jacobian6(m000,m001,m002,m003,m004,m100,m101,m102,m103,m200,m201,m202,m300,m301,m400);
        lam6b = eig(J6);
        lam6br = sort(real(lam6b));
        v6min = min([v6min lam6br(1)]);
        v6max = max([v6max lam6br(6)]);
    else % y-direction
        J6 = jacobian6(m000,m100,m200,m300,m400,m010,m110,m210,m310,m020,m120,m220,m030,m130,m040);
        lam6a = eig(J6);
        lam6ar = sort(real(lam6a));
        v6min = lam6ar(1);
        v6max = lam6ar(6);
        
        J6 = jacobian6(m000,m001,m002,m003,m004,m010,m011,m012,m013,m020,m021,m022,m030,m031,m040);
        lam6b = eig(J6);
        lam6br = sort(real(lam6b));
        v6min = min([v6min lam6br(1)]);
        v6max = max([v6max lam6br(6)]);
    end
end

end
