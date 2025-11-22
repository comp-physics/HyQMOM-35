using Test
using HyQMOM

const HAS_MPI_SIM = try
    using MPI
    true
catch
    false
end

@testset "High-level simulation API" begin
    if !HAS_MPI_SIM
        @test_skip "MPI not available; skipping high-level simulation tests"
    else
        # Ensure MPI is initialized (but do not finalize it here)
        if !MPI.Initialized()
            MPI.Init()
        end

        comm = MPI.COMM_WORLD
        rank = MPI.Comm_rank(comm)
        nprocs = MPI.Comm_size(comm)

        @test nprocs ≥ 1

        # ------------------------------------------------------------------
        # Smoke test for run_simulation (serial or MPI, small problem)
        # ------------------------------------------------------------------
        Nx, Ny, Nz = 6, 4, 1
        tmax = 0.01

        results = HyQMOM.run_simulation(
            Nx = Nx,
            Ny = Ny,
            Nz = Nz,
            tmax = tmax,
            num_workers = nprocs,
            verbose = false,
            enable_plots = false,
            save_output = false,
            homogeneous_z = true,
            debug_output = false,
        )

        @test isa(results, Dict)
        @test haskey(results, :M)
        @test haskey(results, :final_time)
        @test haskey(results, :time_steps)

        if rank == 0
            M = results[:M]
            @test M !== nothing
            @test size(M) == (Nx, Ny, Nz, 35)
            @test results[:final_time] > 0
            @test results[:time_steps] > 0
        else
            @test results[:M] === nothing
        end

        # ------------------------------------------------------------------
        # Lightweight snapshot-mode test for run_simulation_with_snapshots
        # ------------------------------------------------------------------
        Nx2, Ny2, Nz2 = 4, 4, 1
        tmax2 = 0.02
        snapshot_interval = 2
        snapshot_filename = "test_snapshot_small.jld2"

        params = (
            Nx = Nx2,
            Ny = Ny2,
            Nz = Nz2,
            tmax = tmax2,
            Kn = 1.0,
            Ma = 0.0,
            flag2D = 0,
            CFL = 0.5,
            Nmom = 35,
            nnmax = 1_000,
            dtmax = 0.05,
            rhol = 1.0,
            rhor = 0.01,
            T = 1.0,
            r110 = 0.0,
            r101 = 0.0,
            r011 = 0.0,
            symmetry_check_interval = 100,
            homogeneous_z = true,
            enable_memory_tracking = false,
            debug_output = false,
            snapshot_filename = snapshot_filename,
        )

        filename, grid = HyQMOM.run_simulation_with_snapshots(
            params; snapshot_interval = snapshot_interval
        )

        if rank == 0
            @test filename === snapshot_filename
            @test grid !== nothing
            @test isfile(filename)

            using JLD2
            jld = JLD2.jldopen(filename, "r")
            try
                @test haskey(jld, "meta/params")
                @test haskey(jld, "grid")
                @test haskey(jld, "meta/n_snapshots")
                n_snaps = jld["meta/n_snapshots"]
                @test n_snaps ≥ 1
            finally
                close(jld)
                rm(filename; force = true)
            end
        end
    end
end


