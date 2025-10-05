function varargout = mpi_utils(operation, varargin)
%MPI_UTILS Unified MPI utility functions for domain decomposition
%
% Syntax:
%   [Px, Py] = mpi_utils('choose_grid', num_workers)
%   [n_local, i0, i1] = mpi_utils('partition', n, P, r)
%   M_full = mpi_utils('gather_M', M_interior, i0i1, j0j1, Np, Nmom)
%   mpi_utils('send_M', M_interior, i0i1, j0j1, dest_rank)
%
% Operations:
%   'choose_grid' - Choose nearly-square process grid factorization
%   'partition'   - Compute 1D block decomposition for a rank
%   'gather_M'    - Gather moment arrays from all ranks (rank 1 only)
%   'send_M'      - Send moment array to destination rank

    switch operation
        case 'choose_grid'
            [varargout{1}, varargout{2}] = choose_process_grid_impl(varargin{:});
        case 'partition'
            [varargout{1}, varargout{2}, varargout{3}] = block_partition_1d_impl(varargin{:});
        case 'gather_M'
            varargout{1} = gather_M_impl(varargin{:});
        case 'send_M'
            send_M_impl(varargin{:});
        otherwise
            error('mpi_utils:unknownOp', 'Unknown operation: %s', operation);
    end
end

function [Px, Py] = choose_process_grid_impl(nl)
% Choose nearly-square factorization Px*Py=nl where Px and Py are as 
% close as possible to sqrt(nl)

    bestDiff = inf;
    Px = 1; 
    Py = nl;
    
    for p = 1:nl
        if mod(nl, p) == 0
            q = nl / p;
            d = abs(p - q);
            if d < bestDiff
                bestDiff = d;
                Px = p;
                Py = q;
            end
        end
    end
end

function [n_local, i0, i1] = block_partition_1d_impl(n, P, r)
% Compute 1D block decomposition for rank r (0-indexed)
% Divides n points among P processes as evenly as possible
% First (n mod P) ranks get (floor(n/P) + 1) points each

    base = floor(n/P);
    remn = mod(n, P);
    
    if r < remn
        n_local = base + 1;
        i0 = r*(base+1) + 1;
    else
        n_local = base;
        i0 = remn*(base+1) + (r-remn)*base + 1;
    end
    
    i1 = i0 + n_local - 1;
end

function M_full = gather_M_impl(M_interior, i0i1, j0j1, Np, Nmom)
% Gather moment arrays from all ranks to rank 1
% Must be called from rank 1 only, inside spmd block

    % Initialize full array on rank 1
    M_full = zeros(Np, Np, Nmom);
    
    % Place rank 1's data
    M_full(i0i1(1):i0i1(2), j0j1(1):j0j1(2), :) = M_interior;
    
    % Receive from all other ranks
    for src = 2:numlabs
        % Receive data packet: {M_interior, i0i1, j0j1}
        data_packet = labReceive(src);
        blk = data_packet{1};
        i_range = data_packet{2};
        j_range = data_packet{3};
        M_full(i_range(1):i_range(2), j_range(1):j_range(2), :) = blk;
    end
end

function send_M_impl(M_interior, i0i1, j0j1, dest_rank)
% Send moment array to destination rank
% Must be called from non-rank-1 workers, inside spmd block

    % Package data with index information
    data_packet = {M_interior, i0i1, j0j1};
    labSend(data_packet, dest_rank);
end
