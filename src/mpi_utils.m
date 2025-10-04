function varargout = mpi_utils(operation, varargin)
%MPI_UTILS Unified MPI utility functions for domain decomposition
%
% Syntax:
%   [Px, Py] = mpi_utils('choose_grid', num_workers)
%   [n_local, i0, i1] = mpi_utils('partition', n, P, r)
%
% Operations:
%   'choose_grid' - Choose nearly-square process grid factorization
%   'partition'   - Compute 1D block decomposition for a rank

    switch operation
        case 'choose_grid'
            [varargout{1}, varargout{2}] = choose_process_grid_impl(varargin{:});
        case 'partition'
            [varargout{1}, varargout{2}, varargout{3}] = block_partition_1d_impl(varargin{:});
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

