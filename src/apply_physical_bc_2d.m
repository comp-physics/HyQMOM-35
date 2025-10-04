function A = apply_physical_bc_2d(A, decomp, bc)
%APPLY_PHYSICAL_BC_2D Fill halos at global boundaries based on bc.type.
%   A = apply_physical_bc_2d(A, decomp, bc)
%
% Inputs:
%   A      - local array (nx+2h) x (ny+2h) x nv
%   decomp - domain decomposition struct
%   bc     - boundary condition struct with field 'type'
%
% Supported bc.type:
%   'copy' - Neumann-like (copy nearest interior cell)

    h  = decomp.halo;
    nx = decomp.local_size(1);
    ny = decomp.local_size(2);

    if h == 0
        return;
    end
    
    if ~isfield(bc, 'type')
        bc.type = 'copy';
    end

    switch bc.type
        case 'copy'
            % Left boundary (global)
            if decomp.neighbors.left == -1
                for ih = 1:h
                    A(ih, :, :) = A(h+1, :, :);
                end
            end
            % Right boundary (global)
            if decomp.neighbors.right == -1
                for ih = 1:h
                    A(h+nx+ih, :, :) = A(h+nx, :, :);
                end
            end
            % Bottom boundary (global)
            if decomp.neighbors.down == -1
                for ih = 1:h
                    A(:, ih, :) = A(:, h+1, :);
                end
            end
            % Top boundary (global)
            if decomp.neighbors.up == -1
                for ih = 1:h
                    A(:, h+ny+ih, :) = A(:, h+ny, :);
                end
            end
        otherwise
            error('Unknown bc.type: %s', bc.type);
    end
end

