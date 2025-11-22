function Mout = collision35(M,dt,Kn)
% collision35(Min,dt,Kn) elastic BGK collisions
%   Input: 
%       Mi = moments up to 4th order
%       dt = time step
%       Kn = Knudsen number
M000 = M(1);
M100 = M(2);
M200 = M(3);
M010 = M(6);
M020 = M(10);
M001 = M(16);
M002 = M(20);

rho = M000;
umean = M100/rho;
vmean = M010/rho;
wmean = M001/rho;
C200 = M200/rho - umean^2;
C020 = M020/rho - vmean^2;
C002 = M002/rho - wmean^2;

% temperature (2-D)
Theta = (C200 + C020 + C002)/3;
sTheta = sqrt(Theta);

% Maxwellian covariance matrix
CG200 = Theta;
CG020 = Theta;
CG002 = Theta;
CG110 = 0;
CG101 = 0;
CG011 = 0;

% collision time scale
tc = Kn/(rho*sTheta*2); 

% Equilibrium distribution
MG = InitializeM4_35(rho,umean,vmean,wmean,CG200,CG110,CG101,CG020,CG011,CG002);

if size(M) ~= size(MG)
    warning('size mismatch')
end

% BGK model (semi-analytical)
Mout = MG - exp(-dt/tc)*(MG-M);

end