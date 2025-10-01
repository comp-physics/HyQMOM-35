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

    % Left/Right exchange: send/recv halo-width columns over interior rows
    if decomp.neighbors.left ~= -1
        sendbuf = A(h+1:h+h, h+1:h+ny, :); % left interior block width=h
        recvbuf = labSendReceive(decomp.neighbors.left, decomp.neighbors.left, sendbuf);
        A(1:h, h+1:h+ny, :) = recvbuf;
    end
    if decomp.neighbors.right ~= -1
        sendbuf = A(h+nx-h+1:h+nx, h+1:h+ny, :); % right interior block width=h
        recvbuf = labSendReceive(decomp.neighbors.right, decomp.neighbors.right, sendbuf);
        A(h+nx+1:h+nx+h, h+1:h+ny, :) = recvbuf;
    end

    % Up/Down exchange: send/recv halo-width rows over full columns (including updated side halos)
    if decomp.neighbors.down ~= -1
        sendbuf = A(h+1:h+nx, h+1:h+h, :); % bottom interior block height=h
        recvbuf = labSendReceive(decomp.neighbors.down, decomp.neighbors.down, sendbuf);
        A(h+1:h+nx, 1:h, :) = recvbuf;
    end
    if decomp.neighbors.up ~= -1
        sendbuf = A(h+1:h+nx, h+ny-h+1:h+ny, :); % top interior block height=h
        recvbuf = labSendReceive(decomp.neighbors.up, decomp.neighbors.up, sendbuf);
        A(h+1:h+nx, h+ny+1:h+ny+h, :) = recvbuf;
    end
end

