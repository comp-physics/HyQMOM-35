"""
Create MPI Golden Files

This script generates golden files for MPI consistency testing.
It runs simulations with 1, 2, 4, and 8 ranks and saves the results
as reference golden files for regression testing.

Usage:
    mpiexec -n 1 julia --project=. test/create_mpi_goldenfiles.jl
"""

using HyQMOM
using MPI
using Printf

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
nprocs = MPI.Comm_size(comm)

if nprocs != 1
    if rank == 0
        println("ERROR: This script must be run with exactly 1 MPI rank")
        println("Usage: mpiexec -n 1 julia --project=. test/create_mpi_goldenfiles.jl")
    end
    MPI.Finalize()
    exit(1)
end

println("="^70)
println("CREATING MPI GOLDEN FILES")
println("="^70)
println()

# Test configurations
configs = [
    (Np=20, tmax=0.1, name="small"),
    (Np=40, tmax=0.05, name="medium"),
]

golden_dir = joinpath(@__DIR__, "goldenfiles")
mkpath(golden_dir)

for config in configs
    Np = config.Np
    tmax = config.tmax
    name = config.name
    
    println("Configuration: $(name) (Np=$(Np), tmax=$(tmax))")
    println("-"^70)
    
    # Run with 1 rank (reference)
    println("  Running 1-rank simulation...")
    results = run_simulation(
        Np = Np,
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
    filename = joinpath(golden_dir, "mpi_1rank_$(name).bin")
    open(filename, "w") do io
        write(io, Int64(1))        # Number of ranks
        write(io, Int64(Np))       # Grid size
        write(io, Float64(tmax))   # Max time
        write(io, Float64(t_final))# Final time
        write(io, Int64(steps))    # Number of steps
        write(io, M_final)         # Final moments
    end
    
    println("  ✓ Saved: $(filename)")
    println("    Final time: $(t_final), Steps: $(steps)")
    println()
end

println("="^70)
println("GOLDEN FILES CREATED SUCCESSFULLY")
println("="^70)
println()
println("Golden files saved in: $(golden_dir)")
println()
println("To test MPI consistency:")
println("  ./test/test_mpi_goldenfiles.sh")

MPI.Finalize()

