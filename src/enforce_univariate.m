function [S3, S4, H] = enforce_univariate(S3, S4, h2min, s3max)
% enforce_univariate Enforce realizability constraints on univariate moments
%
% Inputs:
%   S3 - third-order standardized moment
%   S4 - fourth-order standardized moment  
%   h2min - minimum allowed value for H = S4 - S3^2 - 1
%   s3max - maximum allowed absolute value for S3
%
% Outputs:
%   S3, S4, H - corrected moments

% Compute and enforce H constraint
H = S4 - S3^2 - 1;
if H <= h2min
    H = h2min;
    S4 = H + S3^2 + 1;
end

% Enforce S3 bounds and update S4 accordingly
if S3 < -s3max || S3 > s3max
    S3 = max(min(S3, s3max), -s3max);
    S4 = H + S3^2 + 1;
end

end
