function Mout = collision15(M,dt,Kn)
% collision12(Min,dt,Kn) elastic BGK collisions
%   Input: 
%       Mi = moments up to 4th order
%       dt = time step
%       Kn = Knudsen number
M00 = M(1);
M10 = M(2);
M20 = M(3);
M01 = M(6);
%M11 = M(7);
M02 = M(10);

rho = M00;
umean = M10/rho;
vmean = M01/rho;
C20 = M20/rho - umean^2;
%C11 = M11/rho - umean*vmean;
C02 = M02/rho - vmean^2;

% temperature (2-D)
Theta = (C20 + C02)/2;
sTheta = sqrt(Theta);

% BGK covariance matrix
%CG20 = (C20+Theta)/2;
%CG11 = C11/2;
%CG02 = (C02+Theta)/2;

% Maxwellian covariance matrix
CG20 = Theta;
CG11 = 0;
CG02 = Theta;

% collision time scale
tc = Kn/(rho*sTheta*2); 

% Equilibrium distribution
MG = InitializeM4_15(rho,umean,vmean,CG20,CG11,CG02); 

% BGK model (semi-analytical)
Mout = MG - exp(-dt/tc)*(MG-M);

end