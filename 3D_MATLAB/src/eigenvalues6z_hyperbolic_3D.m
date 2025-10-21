function [v6min, v6max, Mr] = eigenvalues6z_hyperbolic_3D(M, flag2D, Ma)
% EIGENVALUES6Z_HYPERBOLIC_3D Eigenvalues of 3-D flux Jacobian in z direction
%   [v6min, v6max, Mr] = eigenvalues6z_hyperbolic_3D(M, flag2D, Ma)
%   Inputs:
%       M      - 35-element moment vector
%       flag2D - 2D simulation flag
%       Ma     - Mach number
%   Outputs:
%       v6min, v6max - Min/max eigenvalues for hyperbolicity
%       Mr           - Corrected moment vector
%
% M = [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,...
%      M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,...
%      M031,M012,M112,M013,M022]

Mr = M;

% Extract moments needed for WU and WV planes
m000 = M(1);
m100 = M(2);
m200 = M(3);
m300 = M(4);
m400 = M(5);
m010 = M(6);
m020 = M(10);
m030 = M(13);
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
m021 = M(29);
m031 = M(31);
m012 = M(32);
m013 = M(34);
m022 = M(35);

% WU moments (z-x plane)
J6 = jacobian6(m000, m100, m200, m300, m400, m001, m101, m201, m301, m002, m102, m202, m003, m103, m004);
lam6a = eig(J6);
lam6ar = sort(real(lam6a));
v6min = lam6ar(1);
v6max = lam6ar(6);

% WV moments (z-y plane)
J6 = jacobian6(m000, m010, m020, m030, m040, m001, m011, m012, m013, m002, m021, m022, m003, m013, m004);
lam6b = eig(J6);
lam6br = sort(real(lam6b));
v6min = min([v6min lam6br(1)]);
v6max = max([v6max lam6br(6)]);

% Check for complex eigenvalues in z direction and correct
if max(abs(imag(lam6a))) > 1000*eps || max(abs(imag(lam6b))) > 1000*eps
    % Compute mean velocities
    M000 = M(1);
    M100 = M(2);
    M010 = M(6);
    M001 = M(16);
    umean = M100/M000;
    vmean = M010/M000;
    wmean = M001/M000;
    
    % Compute central and standardized moments
    [C4, S4] = M2CS4_35(M);
    
    % Extract central moments
    C200 = C4(3);
    C300 = C4(4);
    C400 = C4(5);
    C110 = C4(7);
    C210 = C4(8);
    C310 = C4(9);
    C020 = C4(10);
    C120 = C4(11);
    C220 = C4(12);
    C030 = C4(13);
    C130 = C4(14);
    C040 = C4(15);
    C101 = C4(17);
    C201 = C4(18);
    C301 = C4(19);
    C002 = C4(20);
    C102 = C4(21);
    C202 = C4(22);
    C003 = C4(23);
    C103 = C4(24);
    C004 = C4(25);
    C011 = C4(26);
    C111 = C4(27);
    C211 = C4(28);
    C021 = C4(29);
    C121 = C4(30);
    C031 = C4(31);
    C012 = C4(32);
    C112 = C4(33);
    C013 = C4(34);
    C022 = C4(35);
    
    % Extract standardized moments
    S101 = S4(17);
    S300 = S4(4);
    S003 = S4(23);
    S400 = S4(5);
    S202 = S4(22);
    S004 = S4(25);
    S011 = S4(26);
    S030 = S4(13);
    S022 = S4(35);
    S040 = S4(15);
    
    % Force real eigenvalues for WU plane
    if max(abs(imag(lam6a))) > 1000*eps
        S102 = S101*S003;
        S201 = S101*S300;
        S103 = S101*S004;
        S301 = S101*S400;
        s22min = (2 + 5*S003*S101*S300) / 6;
        if S202 < s22min
            S202 = s22min;
        end
        % Central moments from corrected standardized moments
        sC200 = sqrt(max(eps, C200));
        sC002 = sqrt(max(eps, C002));
        
        C201 = S201*sC200^2*sC002;
        C102 = S102*sC200*sC002^2;
        C301 = S301*sC200^3*sC002;
        C103 = S103*sC200*sC002^3;
        C202 = S202*sC200^2*sC002^2;
    end
    
    % Force real eigenvalues for WV plane
    if max(abs(imag(lam6b))) > 1000*eps
        S012 = S011*S003;
        S021 = S011*S030;
        S013 = S011*S004;
        S031 = S011*S040;
        s22min = (2 + 5*S003*S011*S030) / 6;
        if S022 < s22min
            S022 = s22min;
        end
        % Central moments from corrected standardized moments
        sC020 = sqrt(max(eps, C020));
        sC002 = sqrt(max(eps, C002));
        
        C021 = S021*sC020^2*sC002;
        C012 = S012*sC020*sC002^2;
        C031 = S031*sC020^3*sC002;
        C013 = S013*sC020*sC002^3;
        C022 = S022*sC020^2*sC002^2;
    end
    
    % Integer moments from corrected central moments
    M4 = C4toM4_3D(M000, umean, vmean, wmean, C200, C110, C101, C020, C011, C002, C300, ...
                   C210, C201, C120, C111, C102, C030, C021, C012, C003, C400, C310, C301, ...
                   C220, C211, C202, C130, C121, C112, C103, C040, C031, C022, C013, C004);
    
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
    
    % Hyperbolic moments
    Mh = [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,...
          M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,...
          M031,M012,M112,M013,M022];
    
    % Apply realizability constraints
    [~, ~, ~, Mr] = Flux_closure35_and_realizable_3D(Mh, flag2D, Ma);
    
    % Recompute eigenvalues with corrected moments
    m000 = Mr(1);
    m100 = Mr(2);
    m200 = Mr(3);
    m300 = Mr(4);
    m400 = Mr(5);
    m010 = Mr(6);
    m020 = Mr(10);
    m030 = Mr(13);
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
    
    % WU moments
    J6 = jacobian6(m000, m100, m200, m300, m400, m001, m101, m201, m301, m002, m102, m202, m003, m103, m004);
    lam6a = eig(J6);
    lam6ar = sort(real(lam6a));
    v6min = lam6ar(1);
    v6max = lam6ar(6);
    
    % WV moments
    J6 = jacobian6(m000, m010, m020, m030, m040, m001, m011, m012, m013, m002, m021, m022, m003, m013, m004);
    lam6b = eig(J6);
    lam6br = sort(real(lam6b));
    v6min = min([v6min lam6br(1)]);
    v6max = max([v6max lam6br(6)]);
end
end

