function [max_diff, rel_diff, passed] = compare_moments(M1, M2, tolerance, description)
%COMPARE_MOMENTS Compare two moment arrays with tolerance
%   [max_diff, rel_diff, passed] = compare_moments(M1, M2, tolerance, description)
%   
%   Inputs:
%       M1, M2      - Moment arrays to compare
%       tolerance   - Absolute tolerance for comparison
%       description - String describing what is being compared
%   
%   Outputs:
%       max_diff  - Maximum absolute difference
%       rel_diff  - Maximum relative difference
%       passed    - Boolean indicating if comparison passed

if nargin < 4
    description = 'Moments';
end

% Check sizes match
if ~isequal(size(M1), size(M2))
    error('Moment array sizes do not match: %s vs %s', ...
          mat2str(size(M1)), mat2str(size(M2)));
end

% Compute differences
diff_M = abs(M1 - M2);
max_diff = max(diff_M(:));
mean_diff = mean(diff_M(:));

% Compute relative difference (avoid division by zero)
max_val = max(abs(M2(:)));
if max_val > eps
    rel_diff = max_diff / max_val;
else
    rel_diff = max_diff;
end

% Check if passed
passed = max_diff < tolerance;

% Display results if verbose
if nargout == 0 || ~passed
    fprintf('\n%s Comparison:\n', description);
    fprintf('  Max absolute difference: %.6e\n', max_diff);
    fprintf('  Mean absolute difference: %.6e\n', mean_diff);
    fprintf('  Max relative difference: %.6e\n', rel_diff);
    fprintf('  Tolerance: %.6e\n', tolerance);
    if passed
        fprintf('  Status: PASS\n');
    else
        fprintf('  Status: FAIL (exceeds tolerance)\n');
    end
end

end
