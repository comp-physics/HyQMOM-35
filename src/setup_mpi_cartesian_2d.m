function decomp = setup_mpi_cartesian_2d(np, halo)
%SETUP_MPI_CARTESIAN_2D Create a 2D Cartesian domain decomposition over labs
%   decomp = setup_mpi_cartesian_2d(np, halo)
% Inputs:
%   np   - global grid size along x and y (np x np)
%   halo - halo width in cells (typically 1)
% Output struct fields:
%   .np_global   - scalar, global grid size
%   .halo        - scalar, halo width
%   .dims        - [Px, Py] process grid dimensions
%   .coords      - [rx, ry] coordinates of this rank in process grid (0-based)
%   .rank        - this lab index (1-based)
%   .neighbors   - struct with fields: left,right,down,up (lab indices or -1)
%   .local_size  - [nx_local, ny_local] interior size (without halos)
%   .istart_iend - [i0, i1] global inclusive interior index range for x
%   .jstart_jend - [j0, j1] global inclusive interior index range for y
% Notes:
% - Uses an approximately square process grid [Px, Py] with Px*Py=numlabs.
% - Uses block decomposition with remainder cells assigned to lower coords.
    arguments
        np (1,1) {mustBeInteger, mustBePositive}
        halo (1,1) {mustBeInteger, mustBeNonnegative}
    end

    % Discover parallel environment
    nl = numlabs;
    if nl < 1
        error('setup_mpi_cartesian_2d:parallelRequired', 'This function must run inside an spmd block.');
    end

    % Choose process grid dimensions (Px, Py)
    [Px, Py] = mpi_utils('choose_grid', nl);

    % Rank coordinates (0-based)
    r = labindex - 1;
    rx = mod(r, Px);
    ry = floor(r / Px);

    % Compute local interior sizes using block decomposition with remainders
    [nx_local, i0, i1] = mpi_utils('partition', np, Px, rx);
    [ny_local, j0, j1] = mpi_utils('partition', np, Py, ry);

    % Neighbor ranks (1-based), -1 if boundary
    neighbors = struct('left', -1, 'right', -1, 'down', -1, 'up', -1);
    if rx > 0
        neighbors.left = rank_from_coords(rx-1, ry, Px);
    end
    if rx < Px-1
        neighbors.right = rank_from_coords(rx+1, ry, Px);
    end
    if ry > 0
        neighbors.down = rank_from_coords(rx, ry-1, Px);
    end
    if ry < Py-1
        neighbors.up = rank_from_coords(rx, ry+1, Px);
    end

    decomp = struct();
    decomp.np_global = np;
    decomp.halo = halo;
    decomp.dims = [Px, Py];
    decomp.coords = [rx, ry];
    decomp.rank = labindex;
    decomp.neighbors = neighbors;
    decomp.local_size = [nx_local, ny_local];
    decomp.istart_iend = [i0, i1];
    decomp.jstart_jend = [j0, j1];
end

function rank = rank_from_coords(rx, ry, Px)
% Convert (rx,ry) 0-based to 1-based lab index
    rank = ry * Px + rx + 1;
end


