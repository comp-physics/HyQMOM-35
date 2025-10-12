"""
Simpler 3D crossing jets - start with smaller grid and homogeneous_z
"""

using HyQMOM
using MPI

# Initialize MPI
MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

if rank == 0
    println("="^70)
    println("3D Crossing Jets - Simplified Test")
    println("="^70)
end

# Simulation parameters - more conservative
params = (
    Np = 20,              # Smaller grid for testing
    Nz = 10,              # Fewer z-levels
    tmax = 0.02,          # Very short time
    Kn = 1.0,
    Ma = 0.3,             # Lower Mach number for stability
    flag2D = 0,
    CFL = 0.4,            # Conservative CFL
    dx = 1.0/20,
    dy = 1.0/20,
    dz = 1.0/10,
    Nmom = 35,
    nnmax = 50,           # Limit steps
    dtmax = 1e-3,
    
    # Initial condition parameters
    rhol = 1.0,
    rhor = 0.01,
    T = 1.0,
    r110 = 0.0,
    r101 = 0.0,
    r011 = 0.0,
    
    # Diagnostic parameters
    symmetry_check_interval = 100,
    homogeneous_z = true,   # Start with homogeneous (easier)
    debug_output = false,
    enable_memory_tracking = false
)

if rank == 0
    println("\nRunning simplified 3D test...")
    println("  Grid: $(params.Np)×$(params.Np)×$(params.Nz)")
    println("  Ma: $(params.Ma), CFL: $(params.CFL)")
    println("  homogeneous_z: $(params.homogeneous_z)")
    println("  tmax: $(params.tmax)")
end

# Run simulation
M_final, final_time, time_steps, grid = simulation_runner(params)

if rank == 0
    println("\n" * "="^70)
    if isnan(final_time)
        println("❌ Simulation failed with NaN")
        println("This indicates numerical instability")
    else
        println("✅ Simulation Complete!")
        println("="^70)
        println("Final time: $(final_time)")
        println("Time steps: $(time_steps)")
    end
    
    if M_final !== nothing
        println("Result array shape: $(size(M_final))")
        
        # Check for NaNs
        if any(isnan.(M_final))
            println("❌ Result contains NaN values")
        else
            println("✅ No NaN values in result")
            
            # Basic statistics
            rho_min = minimum(M_final[:,:,:,1])
            rho_max = maximum(M_final[:,:,:,1])
            println("Density range: [$(rho_min), $(rho_max)]")
        end
    end
    
    println("="^70)
end

MPI.Finalize()

