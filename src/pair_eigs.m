function [emin, emax] = pair_eigs(JA, JB)
% pair_eigs Compute min/max eigenvalues from two Jacobian matrices
%
% Inputs:
%   JA, JB - Jacobian matrices
%
% Outputs:
%   emin - minimum real eigenvalue from both matrices
%   emax - maximum real eigenvalue from both matrices

% Compute eigenvalues for first matrix
ea = sort(real(eig(JA)));
emin = ea(1);
emax = ea(end);

% Compute eigenvalues for second matrix and update bounds
eb = sort(real(eig(JB)));
emin = min(emin, eb(1));
emax = max(emax, eb(end));

end
