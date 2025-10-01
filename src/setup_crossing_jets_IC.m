function M = setup_crossing_jets_IC(Np, Nmom, rhol, rhor, Ma, T, r110, r101, r011)
% SETUP_CROSSING_JETS_IC Creates initial conditions for crossing jets problem
%
%   M = setup_crossing_jets_IC(Np, Nmom, rhol, rhor, Ma, T, r110, r101, r011)
%
%   Inputs:
%       Np - Grid points per dimension
%       Nmom - Number of moments (35)
%       rhol - High density (jets)
%       rhor - Low density (background)
%       Ma - Mach number
%       T - Temperature
%       r110, r101, r011 - Correlation coefficients
%
%   Output:
%       M - Np x Np x Nmom array of initial moments

% Initialize moment array
M = zeros(Np, Np, Nmom);

% Background parameters (mean velocities at rest)
U0 = 0;
V0 = 0;
W0 = 0;

% Covariance matrix for background
C200 = T;
C020 = T;
C002 = T;
C110 = r110*sqrt(C200*C020);
C101 = r101*sqrt(C200*C002);
C011 = r011*sqrt(C020*C002);

% Initialize background (low density)
Mr = InitializeM4_35(rhor, U0, V0, W0, C200, C110, C101, C020, C011, C002);
for i = 1:Np
    for j = 1:Np
        M(i,j,:) = Mr;
    end
end

% Crossing jets parameters
Uc = Ma / sqrt(2);  % Velocity magnitude for crossing jets
Mt = InitializeM4_35(rhol, -Uc, -Uc, W0, C200, C110, C101, C020, C011, C002);
Mb = InitializeM4_35(rhol,  Uc,  Uc, W0, C200, C110, C101, C020, C011, C002);

% Define jet regions (center 10% of domain on each side)
Csize = floor(0.1*Np);
Mint = Np/2 + 1;
Maxt = Np/2 + 1 + Csize;
Minb = Np/2 - Csize;
Maxb = Np/2;

% Set bottom jet (moving up-right)
for i = Minb:Maxb
    for j = Minb:Maxb
        M(i,j,:) = Mb;
    end
end

% Set top jet (moving down-left)
for i = Mint:Maxt
    for j = Mint:Maxt
        M(i,j,:) = Mt;
    end
end

end

