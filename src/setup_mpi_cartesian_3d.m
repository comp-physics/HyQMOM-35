function decomp = setup_mpi_cartesian_3d(np, nz, halo)
%SETUP_MPI_CARTESIAN_3D Create a 2D Cartesian domain decomposition over labs (no z decomposition)
%   decomp = setup_mpi_cartesian_3d(np, nz, halo)
% Inputs:
%   np   - global grid size along x and y (np x np)
%   nz   - global grid size along z (no decomposition, all ranks have full z)
%   halo - halo width in cells (typically 2)
% Output struct fields:
%   .np_global   - scalar, global grid size in x and y
%   .nz_global   - scalar, global grid size in z
%   .halo        - scalar, halo width
%   .dims        - [Px, Py, 1] process grid dimensions (Pz=1, no z decomposition)
%   .coords      - [rx, ry, 0] coordinates of this rank in process grid (0-based)
%   .rank        - this lab index (1-based)
%   .neighbors   - struct with fields: left,right,down,up (lab indices or -1)
%   .local_size  - [nx_local, ny_local, nz] interior size (without halos)
%   .istart_iend - [i0, i1] global inclusive interior index range for x
%   .jstart_jend - [j0, j1] global inclusive interior index range for y
%   .kstart_kend - [1, nz] global inclusive interior index range for z (always full)
% Notes:
% - Uses an approximately square process grid [Px, Py, 1] with Px*Py=spmdSize.
% - Uses block decomposition in x and y with remainder cells assigned to lower coords.
% - No decomposition in z: all ranks have the full z-dimension (Pz=1)
    arguments
        np (1,1) {mustBeInteger, mustBePositive}
        nz (1,1) {mustBeInteger, mustBePositive}
        halo (1,1) {mustBeInteger, mustBeNonnegative}
    end

    % Discover parallel environment
    nl = spmdSize;
    if nl < 1
        error('setup_mpi_cartesian_3d:parallelRequired', 'This function must run inside an spmd block.');
    end

    % Choose process grid dimensions (Px, Py) - no z decomposition
    [Px, Py] = mpi_utils('choose_grid', nl);
    Pz = 1;  % No decomposition in z

    % Rank coordinates (0-based)
    r = spmdIndex - 1;
    rx = mod(r, Px);
    ry = floor(r / Px);
    rz = 0;  % Always 0 since Pz=1

    % Compute local interior sizes using block decomposition with remainders
    [nx_local, i0, i1] = mpi_utils('partition', np, Px, rx);
    [ny_local, j0, j1] = mpi_utils('partition', np, Py, ry);
    
    % Z dimension: all ranks have the full z extent
    nz_local = nz;
    k0 = 1;
    k1 = nz;

    % Neighbor ranks (1-based), -1 if boundary
    % Only x-y neighbors, no z-neighbors since Pz=1
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
    decomp.nz_global = nz;
    decomp.halo = halo;
    decomp.dims = [Px, Py, Pz];
    decomp.coords = [rx, ry, rz];
    decomp.rank = spmdIndex;
    decomp.neighbors = neighbors;
    decomp.local_size = [nx_local, ny_local, nz_local];
    decomp.istart_iend = [i0, i1];
    decomp.jstart_jend = [j0, j1];
    decomp.kstart_kend = [k0, k1];
end

function rank = rank_from_coords(rx, ry, Px)
% Convert (rx,ry) 0-based to 1-based lab index
    rank = ry * Px + rx + 1;
end


