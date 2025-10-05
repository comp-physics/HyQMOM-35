function M_slice = axis_moment_slice(M, axis)
%AXIS_MOMENT_SLICE Extract 15-element moment vector for jacobian6 computation
%   M_slice = axis_moment_slice(M, axis)
%   Inputs:
%       M    - 35-element moment vector
%       axis - Direction: 1=X (UV plane), 2=Y (VU plane), 3=Z (UW plane), 4=VW plane
%   Outputs:
%       M_slice - 15-element moment vector in jacobian6 order:
%                 [M000, M010, M020, M030, M040, M100, M110, M120, M130,
%                  M200, M210, M220, M300, M310, M400]
%   This eliminates repeated index slicing logic in eigenvalue computations
%   and plotting functions.
%   Examples:
%     % X-direction eigenvalues (UV plane)
%     J6 = jacobian6(axis_moment_slice(M, 1));
%     
%     % Y-direction eigenvalues (VU plane)
%     J6 = jacobian6(axis_moment_slice(M, 2));
%   See also: eigenvalues6_hyperbolic_3D, compute_jacobian_eigenvalues
    % Define moment index mapping for each axis
    switch axis
        case 1  % X direction: UV plane
            % M000, M010, M020, M030, M040, M100, M110, M120, M130, M200, M210, M220, M300, M310, M400
            indices = [1, 6, 10, 13, 15, 2, 7, 11, 14, 3, 8, 12, 4, 9, 5];
            
        case 2  % Y direction: VU plane
            % M000, M100, M200, M300, M400, M010, M110, M210, M310, M020, M120, M220, M030, M130, M040
            indices = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
            
        case 3  % X direction: UW plane
            % M000, M001, M002, M003, M004, M100, M101, M102, M103, M200, M201, M202, M300, M301, M400
            indices = [1, 16, 20, 23, 25, 2, 17, 21, 24, 3, 18, 22, 4, 19, 5];
            
        case 4  % Y direction: VW plane
            % M000, M001, M002, M003, M004, M010, M011, M012, M013, M020, M021, M022, M030, M031, M040
            indices = [1, 16, 20, 23, 25, 6, 26, 32, 34, 10, 29, 35, 13, 31, 15];
            
        otherwise
            error('axis_moment_slice:invalidAxis', ...
                  'Axis must be 1 (X-UV), 2 (Y-VU), 3 (X-UW), or 4 (Y-VW)');
    end
    
    M_slice = M(indices);
end

