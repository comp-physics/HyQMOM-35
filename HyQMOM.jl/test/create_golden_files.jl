"""
Create Golden Files for MPI Testing

This script generates golden reference files for MPI consistency testing.
It runs simulations with 1 rank and saves the results as reference data.

Usage:
    mpiexec -n 1 julia --project=. test/create_golden_files.jl
    mpiexec -n 1 julia --project=. test/create_golden_files.jl small medium
"""

using MPI
using Printf

# Initialize MPI
MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
nprocs = MPI.Comm_size(comm)

if nprocs != 1
    if rank == 0
        println("ERROR: This script must be run with exactly 1 MPI rank")
        println("Usage: mpiexec -n 1 julia --project=. test/create_golden_files.jl")
    end
    MPI.Finalize()
    exit(1)
end

# Load HyQMOM
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using HyQMOM

println("="^70)
println("CREATING MPI GOLDEN FILES")
println("="^70)
println()

# Test configurations
const TEST_CONFIGS = Dict(
    "small" => (Nx=20, Ny=20, tmax=0.1),
    "medium" => (Nx=40, Ny=40, tmax=0.05),
)

# Determine which configs to generate
configs_to_generate = if isempty(ARGS)
    collect(keys(TEST_CONFIGS))
else
    ARGS
end

golden_dir = joinpath(@__DIR__, "goldenfiles")
mkpath(golden_dir)

for config_name in configs_to_generate
    if !haskey(TEST_CONFIGS, config_name)
        println("WARNING: Unknown configuration '$(config_name)', skipping...")
        continue
    end
    
    config = TEST_CONFIGS[config_name]
    Nx = config.Nx
    Ny = config.Ny
    tmax = config.tmax
    
    println("Configuration: $(config_name) (Nx=$(Nx), Ny=$(Ny), tmax=$(tmax))")
    println("-"^70)
    
    # Run with 1 rank (reference)
    println("  Running 1-rank simulation...")
    results = run_simulation(
        Nx = Nx,
        Ny = Ny,
        tmax = tmax,
        num_workers = 1,
        verbose = false,
        Kn = 1.0,
        Ma = 0.0,
        flag2D = 0,
        CFL = 0.5
    )
    
    M_final = results[:M]
    t_final = results[:final_time]
    steps = results[:time_steps]
    
    # Save as golden file
    filename = joinpath(golden_dir, "mpi_1rank_$(config_name).bin")
    open(filename, "w") do io
        write(io, Int64(1))        # Number of ranks
        write(io, Int64(Np))       # Grid size
        write(io, Float64(tmax))   # Max time
        write(io, Float64(t_final))# Final time
        write(io, Int64(steps))    # Number of steps
        write(io, M_final)         # Final moments
    end
    
    println("  ✓ Saved: $(filename)")
    @printf("    Final time: %.6f, Steps: %d\n", t_final, steps)
    println()
end

println("="^70)
println("GOLDEN FILES CREATED SUCCESSFULLY")
println("="^70)
println()
println("Golden files saved in: $(golden_dir)")
println()
println("To test MPI consistency:")
println("  ./test/run_mpi_tests.sh")
println("  or")
println("  mpiexec -n 2 julia --project=. test/test_mpi.jl --golden small")

MPI.Finalize()

