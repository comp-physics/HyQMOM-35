function A = halo_exchange_2d(A, decomp, bc)
%HALO_EXCHANGE_2D Exchange halos for a 2D subdomain with nv variables.
%   A = halo_exchange_2d(A, decomp, bc)
%
% Inputs:
%   A      - local array (nx+2h) x (ny+2h) x nv, interior: A(h+1:h+nx, h+1:h+ny, :)
%   decomp - struct from setup_mpi_cartesian_2d
%   bc     - (optional) boundary condition struct, default: struct('type','copy')
%
% Outputs:
%   A      - updated array with exchanged halos
%
% Notes:
%   - Performs left/right exchange first, then up/down exchange
%   - Applies physical BC at global boundaries before exchange
%   - Corners are filled implicitly after second phase

    h  = decomp.halo;
    nx = decomp.local_size(1);
    ny = decomp.local_size(2);
    
    if nargin < 3
        bc = struct('type', 'copy');
    end

    if h == 0
        return;
    end

    % Apply physical BC to halos at global boundaries before exchange
    A = apply_physical_bc_2d(A, decomp, bc);

    if numlabs == 1
        return;
    end

    % X-direction exchange using non-blocking send/receive
    % Post all sends first, then do receives
    left_neighbor = decomp.neighbors.left;
    right_neighbor = decomp.neighbors.right;
    
    % Send to left neighbor (my leftmost h interior columns)
    if left_neighbor ~= -1
        labSend(A(h+1:h+h, h+1:h+ny, :), left_neighbor);
    end
    
    % Send to right neighbor (my rightmost h interior columns)
    if right_neighbor ~= -1
        labSend(A(h+nx-h+1:h+nx, h+1:h+ny, :), right_neighbor);
    end
    
    % Receive from left neighbor into my left halo
    if left_neighbor ~= -1
        A(1:h, h+1:h+ny, :) = labReceive(left_neighbor);
    end
    
    % Receive from right neighbor into my right halo
    if right_neighbor ~= -1
        A(h+nx+1:h+nx+h, h+1:h+ny, :) = labReceive(right_neighbor);
    end
    
    % Y-direction exchange using non-blocking send/receive
    down_neighbor = decomp.neighbors.down;
    up_neighbor = decomp.neighbors.up;
    
    % Send to down neighbor (my bottommost h interior rows)
    if down_neighbor ~= -1
        labSend(A(h+1:h+nx, h+1:h+h, :), down_neighbor);
    end
    
    % Send to up neighbor (my topmost h interior rows)
    if up_neighbor ~= -1
        labSend(A(h+1:h+nx, h+ny-h+1:h+ny, :), up_neighbor);
    end
    
    % Receive from down neighbor into my bottom halo
    if down_neighbor ~= -1
        A(h+1:h+nx, 1:h, :) = labReceive(down_neighbor);
    end
    
    % Receive from up neighbor into my top halo
    if up_neighbor ~= -1
        A(h+1:h+nx, h+ny+1:h+ny+h, :) = labReceive(up_neighbor);
    end
end

