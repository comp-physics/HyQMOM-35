function varargout = grid_moment_processor(M, func_handle)
%GRID_MOMENT_PROCESSOR Apply any moment function to entire grid
%   This generic utility applies a moment transformation function to every
%   point in a grid, handling the nested loop and memory allocation automatically.
%   Syntax:
%     [C, S] = grid_moment_processor(M, @M2CS4_35)
%     [M5, C5, S5] = grid_moment_processor(M, @Moments5_3D)
%   Inputs:
%     M           - Np x Np x Nmom (2D physical) or Np x Np x Nz x Nmom (3D physical)
%     func_handle - Function handle to apply at each grid point
%                   Must accept a vector of moments and return 1-3 outputs
%   Outputs:
%     varargout   - Variable number of outputs (1-3), each matching input grid dimensions
    
    % Detect dimensionality
    sz = size(M);
    
    if length(sz) == 3
        % 2D physical space: Np x Np x Nmom
        Np = sz(1);
        
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
        
    elseif length(sz) == 4
        % 3D physical space: Np x Np x Nz x Nmom
        Np = sz(1);
        Nz = sz(3);
        
        % Call function once to determine output sizes
        MOM_first = squeeze(M(1,1,1,:));
        [temp_outputs{1:nargout}] = func_handle(MOM_first);
        
        % Pre-allocate output arrays based on first call
        for k = 1:nargout
            out_size = length(temp_outputs{k});
            varargout{k} = zeros(Np, Np, Nz, out_size);
            varargout{k}(1,1,1,:) = temp_outputs{k};
        end
        
        % Process remaining grid points
        for i = 1:Np
            for j = 1:Np
                for kk = 1:Nz
                    if i == 1 && j == 1 && kk == 1
                        continue;  % Already processed
                    end
                    
                    MOM = squeeze(M(i,j,kk,:));
                    [temp_outputs{1:nargout}] = func_handle(MOM);
                    
                    % Store results
                    for k = 1:nargout
                        varargout{k}(i,j,kk,:) = temp_outputs{k};
                    end
                end
            end
        end
        
    else
        error('grid_moment_processor:invalidDims', 'M must be 3D (Np x Np x Nmom) or 4D (Np x Np x Nz x Nmom), got size: [%s]', num2str(sz));
    end
end

