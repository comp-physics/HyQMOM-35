function [n_local, i0, i1] = block_partition_1d(n, P, r)
% BLOCK_PARTITION_1D Compute 1D block decomposition for rank r
%
% Syntax:
%   [n_local, i0, i1] = block_partition_1d(n, P, r)
%
% Description:
%   Divides n points among P processes as evenly as possible, and returns
%   the local size and index range for rank r (0-indexed).
%
%   The first (n mod P) ranks get (floor(n/P) + 1) points each,
%   and the remaining ranks get floor(n/P) points each.
%
% Inputs:
%   n - Total number of points to distribute
%   P - Number of processes
%   r - Rank index (0-based: 0, 1, ..., P-1)
%
% Outputs:
%   n_local - Number of points assigned to rank r
%   i0      - Starting global index for rank r (1-based)
%   i1      - Ending global index for rank r (1-based)
%
% Example:
%   % Distribute 10 points among 3 ranks:
%   [n0, i0_0, i1_0] = block_partition_1d(10, 3, 0)  % n0=4, i0=1, i1=4
%   [n1, i0_1, i1_1] = block_partition_1d(10, 3, 1)  % n1=3, i0=5, i1=7
%   [n2, i0_2, i1_2] = block_partition_1d(10, 3, 2)  % n2=3, i0=8, i1=10

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

