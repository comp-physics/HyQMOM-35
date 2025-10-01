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

    % X-direction exchange: left-right pairs exchange simultaneously
    % Each worker sends right interior edge to right neighbor, receives from left neighbor
    left_neighbor = decomp.neighbors.left;
    right_neighbor = decomp.neighbors.right;
    
    if right_neighbor ~= -1
        % Send my right interior edge to right neighbor
        sendbuf = A(h+nx+1-h:h+nx, h+1:h+ny, :);
        recvbuf = labSendReceive(right_neighbor, right_neighbor, sendbuf);
        % Receive from right neighbor into my right halo
        A(h+nx+1:h+nx+h, h+1:h+ny, :) = recvbuf;
    end
    
    if left_neighbor ~= -1
        % Send my left interior edge to left neighbor
        sendbuf = A(h+1:h+h, h+1:h+ny, :);
        recvbuf = labSendReceive(left_neighbor, left_neighbor, sendbuf);
        % Receive from left neighbor into my left halo
        A(1:h, h+1:h+ny, :) = recvbuf;
    end
    
    % Y-direction exchange: down-up pairs exchange simultaneously
    % Each worker sends top interior edge to up neighbor, receives from down neighbor
    down_neighbor = decomp.neighbors.down;
    up_neighbor = decomp.neighbors.up;
    
    if up_neighbor ~= -1
        % Send my top interior edge to up neighbor
        sendbuf = A(h+1:h+nx, h+ny+1-h:h+ny, :);
        recvbuf = labSendReceive(up_neighbor, up_neighbor, sendbuf);
        % Receive from up neighbor into my top halo
        A(h+1:h+nx, h+ny+1:h+ny+h, :) = recvbuf;
    end
    
    if down_neighbor ~= -1
        % Send my bottom interior edge to down neighbor
        sendbuf = A(h+1:h+nx, h+1:h+h, :);
        recvbuf = labSendReceive(down_neighbor, down_neighbor, sendbuf);
        % Receive from down neighbor into my bottom halo
        A(h+1:h+nx, 1:h, :) = recvbuf;
    end
end

