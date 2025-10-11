function A = halo_exchange_3d(A, decomp, bc)
%HALO_EXCHANGE_3D Exchange halos for a 3D subdomain with nv variables (no z decomposition).
%   A = halo_exchange_3d(A, decomp, bc)
% Inputs:
%   A      - local array (nx+2h) x (ny+2h) x nz x nv, interior: A(h+1:h+nx, h+1:h+ny, :, :)
%   decomp - struct from setup_mpi_cartesian_3d
%   bc     - (optional) boundary condition struct, default: struct('type','copy')
% Outputs:
%   A      - updated array with exchanged halos
% Notes:
%   - Performs left/right exchange first, then up/down exchange
%   - Applies physical BC at global boundaries before exchange
%   - No exchange in z-direction (all ranks have full z extent)
%   - Corners are filled implicitly after second phase
    h  = decomp.halo;
    nx = decomp.local_size(1);
    ny = decomp.local_size(2);
    nz = decomp.local_size(3);
    
    if nargin < 3
        bc = struct('type', 'copy');
    end

    if h == 0
        return;
    end

    % Apply physical BC to halos at global boundaries before exchange
    A = apply_physical_bc_3d(A, decomp, bc);

    if spmdSize == 1
        return;
    end

    % X-direction exchange using non-blocking send/receive
    % Post all sends first, then do receives
    left_neighbor = decomp.neighbors.left;
    right_neighbor = decomp.neighbors.right;
    
    % Send interior boundary data to neighbors (includes all z and all variables)
    if left_neighbor ~= -1
        spmdSend(A(h+1:h+h, h+1:h+ny, :, :), left_neighbor);
    end
    
    if right_neighbor ~= -1
        spmdSend(A(h+nx-h+1:h+nx, h+1:h+ny, :, :), right_neighbor);
    end
    
    % Receive halo data from neighbors
    if left_neighbor ~= -1
        A(1:h, h+1:h+ny, :, :) = spmdReceive(left_neighbor);
    end
    
    if right_neighbor ~= -1
        A(h+nx+1:h+nx+h, h+1:h+ny, :, :) = spmdReceive(right_neighbor);
    end
    
    % Y-direction exchange using non-blocking send/receive
    down_neighbor = decomp.neighbors.down;
    up_neighbor = decomp.neighbors.up;
    
    % Send interior boundary data to neighbors (includes all z and all variables)
    if down_neighbor ~= -1
        spmdSend(A(h+1:h+nx, h+1:h+h, :, :), down_neighbor);
    end
    
    if up_neighbor ~= -1
        spmdSend(A(h+1:h+nx, h+ny-h+1:h+ny, :, :), up_neighbor);
    end
    
    % Receive halo data from neighbors
    if down_neighbor ~= -1
        A(h+1:h+nx, 1:h, :, :) = spmdReceive(down_neighbor);
    end
    
    if up_neighbor ~= -1
        A(h+1:h+nx, h+ny+1:h+ny+h, :, :) = spmdReceive(up_neighbor);
    end
end


