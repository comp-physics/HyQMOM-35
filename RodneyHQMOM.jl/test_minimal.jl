"""
Minimal test to verify the solver works.

This runs a very small simulation (20×20 grid, 2 time steps) to check that
all components are properly connected.
"""

using Printf

# Add src to path
push!(LOAD_PATH, joinpath(@__DIR__, "src"))

# Load MPI first
using MPI
MPI.Init()

# Now load the module
include("src/RodneyHQMOM.jl")
using .RodneyHQMOM

function minimal_test()
    println("="^70)
    println("Minimal Test - Verifying Solver Functionality")
    println("="^70)
    
    # Very small problem for quick test
    Np = 20          # Small grid
    tmax = 0.0001    # Very short time
    nnmax = 2        # Just 2 time steps
    
    # Standard parameters
    Ma = 2.0
    Kn = 0.01
    CFL = 0.5        # Conservative CFL
    flag2D = 0
    
    dx = 1.0 / Np
    dy = 1.0 / Np
    Nmom = 35
    dtmax = 0.001
    
    # IC parameters
    rhol = 1.0
    rhor = 0.5
    T = 1.0
    r110 = 0.0
    r101 = 0.0
    r011 = 0.0
    
    symmetry_check_interval = 1
    
    # Package parameters
    params = (Np=Np, tmax=tmax, Kn=Kn, Ma=Ma, flag2D=flag2D, CFL=CFL,
              dx=dx, dy=dy, Nmom=Nmom, nnmax=nnmax, dtmax=dtmax,
              rhol=rhol, rhor=rhor, T=T, r110=r110, r101=r101, r011=r011,
              symmetry_check_interval=symmetry_check_interval,
              enable_memory_tracking=false)
    
    println("\nTest Configuration:")
    println("  Grid: $(Np)×$(Np)")
    println("  Max steps: $(nnmax)")
    println("  Max time: $(tmax)")
    println()
    
    # Run simulation
    println("Running simulation...")
    start_time = time()
    
    try
        M_final, final_time, time_steps, grid_out = RodneyHQMOM.simulation_runner(params)
        
        elapsed = time() - start_time
        
        comm = MPI.COMM_WORLD
        rank = MPI.Comm_rank(comm)
        
        if rank == 0
            println("\n" * "="^70)
            println("✓ TEST PASSED!")
            println("="^70)
            println("Results:")
            println("  Final time: $(final_time)")
            println("  Time steps: $(time_steps)")
            println("  Wall time: $(elapsed) seconds")
            println("  Final M_final size: $(size(M_final))")
            println("  M_final[1,1,1] (density): $(M_final[1,1,1])")
            println("  M_final range: [$(minimum(M_final)), $(maximum(M_final))]")
            println()
            println("All components working correctly!")
            println("="^70)
            
            return true
        end
        
        return true
        
    catch e
        println("\n" * "="^70)
        println("✗ TEST FAILED!")
        println("="^70)
        println("Error: ", e)
        println()
        println("Stack trace:")
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println()
        end
        println("="^70)
        
        return false
    end
end

# Run test
success = minimal_test()

# Finalize MPI
MPI.Finalize()

# Exit with appropriate code
exit(success ? 0 : 1)
