function [M5, C5, S5] = compute_M5_grid(M)
% compute_M5_grid Compute 5th-order moments on full grid
%
% Input:
%   M - moment array (Np x Np x Nmom)
%
% Outputs:
%   M5 - 5th-order moments (Np x Np x Nmom5)
%   C5 - 5th-order central moments (Np x Np x Nmom5)
%   S5 - 5th-order standardized moments (Np x Np x Nmom5)

[Np, ~, ~] = size(M);
Nmom5 = 21; % Fixed for this application

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
