function [v6min,v6max,Mr] = eigenvalues6x_hyperbolic_3D(M,flag2D,Ma)
% eigenvalues of 3-D flux Jacobian in x direction
%
% M = [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,...
%      M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,...
%      M031,M012,M112,M013,M022]
%
Mr = M;
%
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
% UV  moments
J6 = jacobian6(m000,m010,m020,m030,m040,m100,m110,m120,m130,m200,m210,m220,m300,m310,m400);
lam6a = eig(J6);
lam6ar = sort(real(lam6a));
v6min = lam6ar(1);
v6max = lam6ar(6);
% UW  moments
J6 = jacobian6(m000,m001,m002,m003,m004,m100,m101,m102,m103,m200,m201,m202,m300,m301,m400);
lam6b = eig(J6);
lam6br = sort(real(lam6b));
v6min = min([v6min lam6br(1)]);
v6max = max([v6max lam6br(6)]);
%
% return
%
% check for complex eigenvalues in x direction and correct
%if lam6a(1) == lam6a(2) || lam6a(5) == lam6a(6) || lam6b(1) == lam6b(2) || lam6b(5) == lam6b(6)
if max(abs(imag(lam6a))) > 1000*eps || max(abs(imag(lam6b))) > 1000*eps % any(imag(lam6a),'all') || any(imag(lam6b),'all') 
    %  warning('complex eigenvalues')
    %  eig(J6)
    % compute mean velocities
    M000 = M(1);
    M100 = M(2);
    M010 = M(6);
    M001 = M(16);
    umean = M100/M000;
    vmean = M010/M000;
    wmean = M001/M000; 
    %
    % compute central and standardized moments and closures from 35 known moments
    [C4,S4] = M2CS4_35(M);
    %
    C200=C4(3);
    C300=C4(4);
    C400=C4(5);
    C110=C4(7);
    %
    C210=C4(8);
    C310=C4(9);
    C020=C4(10);
    C120=C4(11);
    C220=C4(12);
    C030=C4(13);
    C130=C4(14);
    C040=C4(15);
    %
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
    %
    C031=C4(31);
    C012=C4(32);
    C112=C4(33);
    C013=C4(34);
    C022=C4(35);
    %
    S110 = S4(7);
    S300 = S4(4);
    S030 = S4(13);
    S400 = S4(5);
    S220 = S4(12);
    S040 = S4(15);
    S101 = S4(17);
    S202 = S4(22);
    S003 = S4(23);
    S004 = S4(25);
    % force real eigenvalues
    if max(abs(imag(lam6a))) > 1000*eps %lam6a(1) == lam6a(2) || lam6a(5) == lam6a(6)
        S120 = S110*S030;
        S210 = S110*S300;
        S130 = S110*S040;
        S310 = S110*S400;
        %s22min = max((2+6*S030*S110*S300-S300^2),(2+6*S030*S110*S300-S030^2))/6;
        s22min = (2+5*S030*S110*S300)/6;
        if S220 < s22min
            S220 = s22min;
        end
        % central moments from corrected standized moments
        sC200 = sqrt(max(eps,C200));
        sC020 = sqrt(max(eps,C020));
        %
        C210 = S210*sC200^2*sC020;
        C120 = S120*sC200*sC020^2;
        C310 = S310*sC200^3*sC020;
        C130 = S130*sC200*sC020^3;
        C220 = S220*sC200^2*sC020^2;
    end
    if max(abs(imag(lam6b))) > 1000*eps %lam6b(1) == lam6b(2) || lam6b(5) == lam6b(6)
        S102 = S101*S003;
        S201 = S101*S300;
        S103 = S101*S004;
        S301 = S101*S400;
        %s22min = max((2+6*S030*S110*S300-S300^2),(2+6*S030*S110*S300-S030^2))/6;
        s22min = (2+5*S003*S101*S300)/6;
        if S202 < s22min
            S202 = s22min;
        end
        % central moments from corrected standized moments
        sC200 = sqrt(max(eps,C200));
        sC002 = sqrt(max(eps,C002));
        %
        C201 = S201*sC200^2*sC002;
        C102 = S102*sC200*sC002^2;
        C301 = S301*sC200^3*sC002;
        C103 = S103*sC200*sC002^3;
        C202 = S202*sC200^2*sC002^2;
    end
    %
    % integer moments from corrected central moments
    M4 = C4toM4_3D(M000,umean,vmean,wmean,C200,C110,C101,C020,C011,C002,C300,...
                   C210,C201,C120,C111,C102,C030,C021,C012,C003,C400,C310,C301,...
                   C220,C211,C202,C130,C121,C112,C103,C040,C031,C022,C013,C004);
    %
    M000 = M4(1,1,1);
    M100 = M4(2,1,1);
    M010 = M4(1,2,1);
    M001 = M4(1,1,2);
    %
    M200 = M4(3,1,1);
    M110 = M4(2,2,1);
    M101 = M4(2,1,2);
    M020 = M4(1,3,1);
    M011 = M4(1,2,2);
    M002 = M4(1,1,3);
    %
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
    %
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
    % hyperbolic moments
    Mh = [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,...
          M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,...
          M031,M012,M112,M013,M022];
    % realizable moments
    [~,~,~,Mr] = Flux_closure35_and_realizable_3D(Mh,flag2D,Ma);
    %
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
    % UV  moments
    J6 = jacobian6(m000,m010,m020,m030,m040,m100,m110,m120,m130,m200,m210,m220,m300,m310,m400);
    lam6a = eig(J6);
    lam6ar = sort(real(lam6a));
    v6min = lam6ar(1);
    v6max = lam6ar(6);
    % UW  moments
    J6 = jacobian6(m000,m001,m002,m003,m004,m100,m101,m102,m103,m200,m201,m202,m300,m301,m400);
    lam6b = eig(J6);
    lam6br = sort(real(lam6b));
    v6min = min([v6min lam6br(1)]);
    v6max = max([v6max lam6br(6)]);
end
end