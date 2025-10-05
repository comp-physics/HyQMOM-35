function varargout = grid_moment_processor(M, func_handle)
%GRID_MOMENT_PROCESSOR Apply any moment function to entire grid
%
%   This generic utility applies a moment transformation function to every
%   point in a grid, handling the nested loop and memory allocation automatically.
%
%   Syntax:
%     [C, S] = grid_moment_processor(M, @M2CS4_35)
%     [M5, C5, S5] = grid_moment_processor(M, @Moments5_3D)
%
%   Inputs:
%     M           - Np x Np x Nmom array of raw moments
%     func_handle - Function handle to apply at each grid point
%                   Must accept a vector of moments and return 1-3 outputs
%
%   Outputs:
%     varargout   - Variable number of outputs (1-3), each as Np x Np x N array

    [Np, ~, ~] = size(M);
    
    % Call function once to determine output sizes
    MOM_first = squeeze(M(1,1,:));
    [temp_outputs{1:nargout}] = func_handle(MOM_first);
    
    % Pre-allocate output arrays based on first call
    for k = 1:nargout
        out_size = length(temp_outputs{k});
        varargout{k} = zeros(Np, Np, out_size);
        varargout{k}(1,1,:) = temp_outputs{k};
    end
    
    % Process remaining grid points
    for i = 1:Np
        for j = 1:Np
            if i == 1 && j == 1
                continue;  % Already processed
            end
            
            MOM = squeeze(M(i,j,:));
            [temp_outputs{1:nargout}] = func_handle(MOM);
            
            % Store results
            for k = 1:nargout
                varargout{k}(i,j,:) = temp_outputs{k};
            end
        end
    end
end

