function [realizable, violations] = verify_realizability(M, tolerance)
%VERIFY_REALIZABILITY Check realizability constraints for moment vector
%   [realizable, violations] = verify_realizability(M, tolerance)
%   
%   Inputs:
%       M         - Moment vector (35 x 1) or array (Nx x Ny x 35)
%       tolerance - Tolerance for constraint violations (default: 1e-6)
%   
%   Outputs:
%       realizable - Boolean indicating if moments are realizable
%       violations - Struct with information about any violations

if nargin < 2
    tolerance = 1e-6;
end

violations = struct();
violations.count = 0;
violations.details = {};

% Handle grid arrays
if ndims(M) == 3
    [Nx, Ny, Nmom] = size(M);
    all_realizable = true;
    
    for i = 1:Nx
        for j = 1:Ny
            M_vec = squeeze(M(i, j, :));
            [real_ij, viol_ij] = check_single_moment(M_vec, tolerance);
            if ~real_ij
                all_realizable = false;
                violations.count = violations.count + viol_ij.count;
                for k = 1:length(viol_ij.details)
                    violations.details{end+1} = sprintf('[%d,%d] %s', ...
                        i, j, viol_ij.details{k});
                end
            end
        end
    end
    realizable = all_realizable;
else
    [realizable, violations] = check_single_moment(M, tolerance);
end

end

function [realizable, violations] = check_single_moment(M, tolerance)
%CHECK_SINGLE_MOMENT Check realizability for a single moment vector

violations = struct();
violations.count = 0;
violations.details = {};

% Basic positivity: M000 > 0
if M(1) <= 0
    violations.count = violations.count + 1;
    violations.details{end+1} = sprintf('M000 = %.6e <= 0', M(1));
end

% Compute standardized moments
[C4, S4] = M2CS4_35(M);

% Check variances are positive
variances = [C4(3), C4(10), C4(20)];  % C200, C020, C002
var_names = {'C200', 'C020', 'C002'};
for i = 1:3
    if variances(i) < -tolerance
        violations.count = violations.count + 1;
        violations.details{end+1} = sprintf('%s = %.6e < 0', ...
            var_names{i}, variances(i));
    end
end

% Check S2 realizability: S2 = 1 + 2*S110*S101*S011 - (S110^2 + S101^2 + S011^2) >= 0
S110 = S4(7);
S101 = S4(17);
S011 = S4(26);
S2 = 1 + 2*S110*S101*S011 - (S110^2 + S101^2 + S011^2);

if S2 < -tolerance
    violations.count = violations.count + 1;
    violations.details{end+1} = sprintf('S2 = %.6e < 0 (realizability violated)', S2);
end

% Check correlation bounds: |S110|, |S101|, |S011| <= 1
if abs(S110) > 1 + tolerance
    violations.count = violations.count + 1;
    violations.details{end+1} = sprintf('|S110| = %.6e > 1', abs(S110));
end
if abs(S101) > 1 + tolerance
    violations.count = violations.count + 1;
    violations.details{end+1} = sprintf('|S101| = %.6e > 1', abs(S101));
end
if abs(S011) > 1 + tolerance
    violations.count = violations.count + 1;
    violations.details{end+1} = sprintf('|S011| = %.6e > 1', abs(S011));
end

realizable = (violations.count == 0);

end
