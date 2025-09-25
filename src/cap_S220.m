function Sr = cap_S220(Sii, Sij, Hii, Hjj, S3i, S3j)
% cap_S220 Apply S220-type realizability capping
%
% Inputs:
%   Sii - second-order moment in first direction
%   Sij - cross moment to be capped
%   Hii, Hjj - H values for both directions  
%   S3i, S3j - third-order moments for both directions
%
% Output:
%   Sr - capped cross moment

A = sqrt((Hii + S3i^2) * (Hjj + S3j^2));
Sr = realizability_engine('S220', Sii, Sij, A);

end
