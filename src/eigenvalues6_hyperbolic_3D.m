function [v6min,v6max,Mr] = eigenvalues6_hyperbolic_3D(M,direction,flag2D,Ma)
% eigenvalues6_hyperbolic_3D - Unified eigenvalue computation for 3D flux Jacobian
%
% Inputs:
%   M         - 35-component moment vector
%   direction - 'x' or 'y' for direction of flux Jacobian
%   flag2D    - flag for 2D case
%   Ma        - Mach number
%
% Outputs:
%   v6min     - minimum eigenvalue
%   v6max     - maximum eigenvalue
%   Mr        - corrected moment vector (if realizability issues)
%
% M = [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,...
%      M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,...
%      M031,M012,M112,M013,M022]
%
Mr = M;
idx = moment_indices();

% Compute eigenvalues using helper
[v6min, v6max, lam6a, lam6b] = compute_eigenvalues_from_moments(M, direction);

% Check for complex eigenvalues and correct if needed
if max(abs(imag(lam6a))) > 1000*eps || max(abs(imag(lam6b))) > 1000*eps
    % Get mean velocities
    [umean, vmean, wmean] = means_from_M(M);
    M000 = M(idx.m000);
    
    % Get central and standardized moments
    [C4,S4] = M2CS4_35(M);
    [sC200, sC020, sC002] = sC_from_C4(C4);
    
    % Extract needed standardized moments
    S110 = S4(7); S101 = S4(17); S011 = S4(26);
    S300 = S4(idx.S300); S030 = S4(idx.S030); S003 = S4(idx.S003);
    S400 = S4(idx.S400); S040 = S4(idx.S040); S004 = S4(idx.S004);
    S220 = S4(12); S202 = S4(22); S022 = S4(35);
    
    % Store needed central moments
    C200 = C4(idx.C200); C020 = C4(idx.C020); C002 = C4(idx.C002);
    C110 = C4(7); C101 = C4(17); C011 = C4(26);
    C220 = C4(12); C202 = C4(22); C022 = C4(35);
    
    % Direction-specific corrections using helper function
    need_correction = false;
    
    if strcmpi(direction, 'x')
        if max(abs(imag(lam6a))) > 1000*eps
            C220 = correct_cross_moment(S400, S040, S300, S030, S110, S220, C200*C020, C110^2, sC200*sC020);
            need_correction = true;
        end
        if max(abs(imag(lam6b))) > 1000*eps
            C202 = correct_cross_moment(S400, S004, S300, S003, S101, S202, C200*C002, C101^2, sC200*sC002);
            need_correction = true;
        end
    else % y-direction
        if max(abs(imag(lam6a))) > 1000*eps
            C220 = correct_cross_moment(S400, S040, S300, S030, S110, S220, C200*C020, C110^2, sC200*sC020);
            need_correction = true;
        end
        if max(abs(imag(lam6b))) > 1000*eps
            C022 = correct_cross_moment(S040, S004, S030, S003, S011, S022, C020*C002, C011^2, sC020*sC002);
            need_correction = true;
        end
    end
    
    if need_correction
        % Reconstruct moments with corrections
        M4 = C4toM4_3D(M000,umean,vmean,wmean,C200,C110,C101,C020,C011,C002,...
            C4(4),C4(8),C4(18),C4(11),C4(27),C4(21),C4(13),C4(29),C4(32),C4(23),...
            C4(5),C4(9),C4(19),C220,C4(28),C202,C4(14),C4(30),C4(33),C4(24),...
            C4(15),C4(31),C022,C4(34),C4(25));
        
        % Pack corrected moments vector
        Mh = zeros(35,1);
        for i = 1:5
            for j = 1:5
                for k = 1:5
                    if i+j+k <= 6 && i+j+k >= 2
                        % Map 3D indices to linear index in 35-element vector
                        Mh(moment_field_order(i,j,k)) = M4(i,j,k);
                    end
                end
            end
        end
        
        % Apply realizability corrections
        [~,~,~,Mr] = Flux_closure35_and_realizable_3D(Mh,flag2D,Ma);
        
        % Recompute eigenvalues with corrected moments
        [v6min, v6max] = compute_eigenvalues_from_moments(Mr, direction);
    end
end

end
