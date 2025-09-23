function [v6min,v6max,Mr] = eigenvalues6y_hyperbolic(M,flagreal,Ma)
% eigenvalues of flux Jacobian in y direction
%
% M = [M00,M10,M20,M30,M40,M01,M11,M21,M31,M02,M12,M22,M03,M13,M04]
%
Mr = M;
%
m00 = M(1);
m10 = M(2);
m20 = M(3);
m30 = M(4);
m40 = M(5);
m01 = M(6);
m11 = M(7);
m21 = M(8);
m31 = M(9);
m02 = M(10);
m12 = M(11);
m22 = M(12);
m03 = M(13);
m13 = M(14);
m04 = M(15);
J6 = jacobian6(m00,m10,m20,m30,m40,m01,m11,m21,m31,m02,m12,m22,m03,m13,m04);
lam6 = eig(J6);
lam6 = sort(real(lam6));
v6min = lam6(1);
v6max = lam6(6);
% check for complex eigenvalues in y direction
if lam6(1) == lam6(2) ||  lam6(5) == lam6(6)
    % compute mean velocities
    M00 = M(1);
    M10 = M(2);
    M01 = M(6);
    umean = M10/M00;
    vmean = M01/M00;
    %
    % compute central and standardized moments and closures from 15 known moments
    [C4,S4] = M2CS4_15(M);
    %
    C20 = C4(3);
    C02 = C4(10);
    %
    S11 = S4(7);
    S30 = S4(4);
    % S21 = S4(8);
    % S12 = S4(11);
    S03 = S4(13);
    S40 = S4(5);
    % S31 = S4(9);
    S22 = S4(12);
    % S13 = S4(14);
    S04 = S4(15);
    % force real eigenvalues
    S12 = S11*S03;
    S21 = S11*S30;
    S13 = S11*S04;
    S31 = S11*S40;
    %s22min = max((2+6*S03*S11*S30-S30^2),(2+6*S03*S11*S30-S03^2))/6;
    s22min = (2+5*S03*S11*S30)/6;
    if S22 < s22min
        S22 = s22min;
    end
    % central moments from standized moments
    sC20 = sqrt(C20);
    sC02 = sqrt(C02);
    C11 = S11*sC20*sC02;
    C30 = S30*sC20^3; 
    C21 = S21*sC20^2*sC02;
    C12 = S12*sC20*sC02^2;
    C03 = S03*sC02^3;
    C40 = S40*sC20^4; 
    C31 = S31*sC20^3*sC02; 
    C22 = S22*sC20^2*sC02^2;
    C13 = S13*sC20*sC02^3;
    C04 = S04*sC02^4;
    %
    % integer moments from central moments
    M4 = C4toM4(M00,umean,vmean,C20,C11,C02,C30,C21,C12,C03,C40,C31,C22,C13,C04);
    %
    M00 = M4(1,1);
    M10 = M4(2,1);
    M01 = M4(1,2);
    M20 = M4(3,1);
    M11 = M4(2,2);
    M02 = M4(1,3);
    M30 = M4(4,1);
    M21 = M4(3,2);
    M12 = M4(2,3);
    M03 = M4(1,4);
    M40 = M4(5,1); 
    M31 = M4(4,2);
    M22 = M4(3,3);
    M13 = M4(2,4);
    M04 = M4(1,5);
    % hyperbolic and realizable moments
    Mh = [M00,M10,M20,M30,M40,M01,M11,M21,M31,M02,M12,M22,M03,M13,M04];
    [~,~,Mr] = Flux_closure15_and_realizable(Mh,flagreal,Ma);
    %
    m00 = Mr(1);
    m10 = Mr(2);
    m20 = Mr(3);
    m30 = Mr(4);
    m40 = Mr(5);
    m01 = Mr(6);
    m11 = Mr(7);
    m21 = Mr(8);
    m31 = Mr(9);
    m02 = Mr(10);
    m12 = Mr(11);
    m22 = Mr(12);
    m03 = Mr(13);
    m13 = Mr(14);
    m04 = Mr(15);
    J6 = jacobian6(m00,m10,m20,m30,m40,m01,m11,m21,m31,m02,m12,m22,m03,m13,m04);
    lam6 = eig(J6);
    lam6 = sort(real(lam6));
    v6min = lam6(1);
    v6max = lam6(6);
end
end