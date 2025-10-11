function [S3, S4, H] = enforce_univariate(S3, S4, h2min, s3max)
% enforce_univariate enforces realizability bounds on univariate moments
% Ensures H = S4 - S3^2 - 1 >= h2min and |S3| <= s3max
% Inputs:
%   S3 - third standardized moment
%   S4 - fourth standardized moment
%   h2min - minimum allowed H value
%   s3max - maximum allowed |S3| value
% Outputs:
%   S3 - corrected third moment
%   S4 - corrected fourth moment  
%   H - corrected H value
% Compute H
H = S4 - S3^2 - 1;

% Enforce minimum H
if H <= h2min
    H = h2min;
    S4 = H + S3^2 + 1;
end

% Enforce S3 bounds
if S3 < -s3max
    S3 = -s3max;
    S4 = H + S3^2 + 1;
elseif S3 > s3max
    S3 = s3max;
    S4 = H + S3^2 + 1;
end

end

