#!/usr/bin/env julia
"""
Simple example: Run a small HyQMOM simulation

This is a minimal example showing how to run a HyQMOM simulation
programmatically using the low-level API (simulation_runner).

For most use cases, prefer the high-level run_simulation() API
shown in run_with_plots.jl.

# Usage
```bash
cd HyQMOM.jl
julia --project=. examples/run_simple.jl
```

# With MPI
```bash
cd HyQMOM.jl
mpiexec -n 2 julia --project=. examples/run_simple.jl
```
"""

using MPI
using HyQMOM
using Printf

# Initialize MPI (required for simulation_runner)
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
nprocs = MPI.Comm_size(comm)

if rank == 0
    println("="^70)
    println("HyQMOM Simple Example - Low-Level API")
    println("="^70)
    println("\nThis example uses simulation_runner() directly.")
    println("For a simpler API, see run_with_plots.jl\n")
end

# Set up parameters for a small test case
# Note: simulation_runner requires all parameters as a named tuple
params = (
    # Grid parameters
    Np = 40,              # Global grid size: 40x40
    dx = 1.0/40,          # Grid spacing (domain is [-0.5, 0.5])
    dy = 1.0/40,
    
    # Time parameters
    tmax = 0.1,           # Maximum simulation time
    dtmax = 1.0,          # Maximum time step (actual dt computed adaptively)
    nnmax = 20000000,     # Maximum number of time steps
    CFL = 0.5,            # CFL number for stability
    
    # Physical parameters
    Kn = 1.0,             # Knudsen number (ratio of mean free path to length scale)
    Ma = 0.0,             # Mach number (ratio of flow speed to sound speed)
    flag2D = 0,           # 0 = 3D simulation, 1 = 2D simulation
    
    # Moment system
    Nmom = 35,            # Number of moments (fixed for 3D HyQMOM)
    
    # Initial condition parameters (crossing jets)
    rhol = 1.0,           # Density in jets (high)
    rhor = 0.01,          # Density in background (low)
    T = 1.0,              # Temperature
    r110 = 0.0,           # Velocity correlations
    r101 = 0.0,
    r011 = 0.0,
    
    # Diagnostic parameters
    symmetry_check_interval = 10,
    enable_memory_tracking = false,
    debug_output = false
)

if rank == 0
    println("Simulation parameters:")
    println("  Grid size: $(params.Np) x $(params.Np)")
    println("  MPI ranks: $nprocs")
    println("  Final time: $(params.tmax)")
    println("  Knudsen number: $(params.Kn)")
    println("  Mach number: $(params.Ma)")
    println("  CFL number: $(params.CFL)")
    println()
    println("Running simulation...")
    flush(stdout)
end

# Run simulation using low-level API
# Returns: M_final, final_time, time_steps, grid
# Note: M_final is only non-nothing on rank 0
M_final, t_final, steps, grid = HyQMOM.simulation_runner(params)

# Display results (only on rank 0)
if rank == 0
    println("\n" * "="^70)
    println("Simulation complete!")
    println("="^70)
    println("  Final time: $t_final")
    println("  Time steps: $steps")
    
    if !isnothing(M_final)
        println("  Result size: $(size(M_final))")
        
        # Extract density (first moment)
        ρ = M_final[:, :, 1]
        @printf("  Density range: [%.6e, %.6e]\n", minimum(ρ), maximum(ρ))
        
        # Check for numerical issues
        if any(isnan, M_final)
            println("  ⚠ WARNING: NaN values detected!")
        elseif any(isinf, M_final)
            println("  ⚠ WARNING: Inf values detected!")
        elseif any(ρ .<= 0)
            println("  ⚠ WARNING: Non-positive densities detected!")
        else
            println("  ✓ All moments are finite and density is positive")
        end
        
        # Basic conservation check
        # Initial mass: integrate ρ over domain
        # For crossing jets: roughly (Np^2 * (rhol + rhor) / 2) * dx * dy
        dx = grid.dx
        dy = grid.dy
        initial_mass_estimate = params.Np * params.Np * (params.rhol + params.rhor) / 2 * dx * dy
        final_mass = sum(ρ) * dx * dy
        mass_error = abs(final_mass - initial_mass_estimate) / initial_mass_estimate
        
        @printf("\n  Mass Conservation:\n")
        @printf("    Initial (estimate): %.6e\n", initial_mass_estimate)
        @printf("    Final:              %.6e\n", final_mass)
        @printf("    Relative error:     %.2e\n", mass_error)
        
        if mass_error < 1e-10
            println("    ✓ Mass is well conserved!")
        elseif mass_error < 1e-6
            println("    ✓ Mass is reasonably conserved")
        else
            println("    ⚠ WARNING: Large mass conservation error")
        end
    end
    
    println("="^70)
end

# Finalize MPI
MPI.Finalize()

