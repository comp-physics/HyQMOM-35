function [C, S] = compute_CS_grid(M)
% compute_CS_grid computes central and standardized moments for entire grid
% 
% Inputs:
%   M - Np x Np x Nmom array of raw moments
%
% Outputs:
%   C - Np x Np x Nmom array of central moments
%   S - Np x Np x Nmom array of standardized moments

[Np, ~, Nmom] = size(M);
C = zeros(Np, Np, Nmom);
S = zeros(Np, Np, Nmom);

for i = 1:Np
    for j = 1:Np
        MOM = squeeze(M(i,j,:));
        [CC, SS] = M2CS4_35(MOM);
        C(i,j,:) = CC;
        S(i,j,:) = SS;
    end
end

end

