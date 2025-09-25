function [C, S] = compute_CS_grid(M)
% compute_CS_grid Compute central and standardized moments on full grid
%
% Input:
%   M - moment array (Np x Np x Nmom)
%
% Outputs:
%   C - central moments (Np x Np x Nmom)
%   S - standardized moments (Np x Np x Nmom)

[Np, ~, Nmom] = size(M);
C = zeros(Np, Np, Nmom);
S = zeros(Np, Np, Nmom);

parfor i = 1:Np
    for j = 1:Np
        MOM = squeeze(M(i,j,:));
        [CC, SS] = M2CS4_35(MOM);
        C(i,j,:) = CC;
        S(i,j,:) = SS;
    end
end

end
