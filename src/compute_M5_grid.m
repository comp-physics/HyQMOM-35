function [M5, C5, S5] = compute_M5_grid(M)
% compute_M5_grid computes 5th-order moments for entire grid
%
% Inputs:
%   M - Np x Np x Nmom array of raw moments (35 moments)
%
% Outputs:
%   M5 - Np x Np x Nmom5 array of 5th-order raw moments (21 moments)
%   C5 - Np x Np x Nmom5 array of 5th-order central moments
%   S5 - Np x Np x Nmom5 array of 5th-order standardized moments

[Np, ~, Nmom] = size(M);
Nmom5 = 21;  % Number of 5th-order moments
M5 = zeros(Np, Np, Nmom5);
C5 = zeros(Np, Np, Nmom5);
S5 = zeros(Np, Np, Nmom5);

for i = 1:Np
    for j = 1:Np
        MOM = squeeze(M(i,j,:));
        [MM5, CC5, SS5] = Moments5_3D(MOM);
        M5(i,j,:) = MM5;
        C5(i,j,:) = CC5;
        S5(i,j,:) = SS5;
    end
end

end

