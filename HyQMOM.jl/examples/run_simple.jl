#!/usr/bin/env julia
"""
Simple example: Run a small HyQMOM simulation

This is a minimal example showing how to run a HyQMOM simulation
programmatically (without command-line arguments).

# Usage
```bash
cd HyQMOM.jl/examples
julia --project=.. run_simple.jl
```

# With MPI
```bash
mpiexec -n 2 julia --project=.. run_simple.jl
```
"""

using HyQMOM
using Printf

println("="^70)
println("HyQMOM Simple Example")
println("="^70)

# Set up parameters for a small test case
params = (
    Np = 40,              # Grid size: 40x40
    tmax = 0.1,           # Final time
    Kn = 1.0,             # Knudsen number
    Ma = 0.0,             # Mach number
    flag2D = 0,           # 3D simulation
    CFL = 0.5,            # CFL number
    dx = 1.0/40,          # Grid spacing
    dy = 1.0/40,
    Nmom = 35,            # Number of moments
    nnmax = 20000000,     # Max time steps
    dtmax = 1.0,          # Max time step
    rhol = 1.0,           # Left density
    rhor = 0.01,          # Right density
    T = 1.0,              # Temperature
    r110 = 0.0,           # Correlations
    r101 = 0.0,
    r011 = 0.0,
    symmetry_check_interval = 10,
    enable_memory_tracking = false
)

println("\nSimulation parameters:")
println("  Grid size: $(params.Np) x $(params.Np)")
println("  Final time: $(params.tmax)")
println("  Knudsen number: $(params.Kn)")
println("  Mach number: $(params.Ma)")
println("  CFL number: $(params.CFL)")
println()

# Run simulation
println("Running simulation...")
M_final, t_final, steps, grid = HyQMOM.simulation_runner(params)

# Display results
println("\nSimulation complete!")
println("  Final time: $t_final")
println("  Time steps: $steps")

if !isnothing(M_final)
    println("  Final moment field size: $(size(M_final))")
    println("  Density range: [$(minimum(M_final[:,:,1])), $(maximum(M_final[:,:,1]))]")
    
    # Basic conservation check
    initial_mass = params.Np * params.Np * (params.rhol + params.rhor) / 2
    final_mass = sum(M_final[:,:,1]) * params.dx * params.dy
    mass_error = abs(final_mass - initial_mass) / initial_mass
    
    @printf("  Mass conservation error: %.2e\n", mass_error)
    
    if mass_error < 1e-10
        println("  OK Mass is conserved!")
    else
        println("  WARNING: Mass conservation error is large")
    end
else
    println("  (Results on rank 0 only)")
end

println("\n" * "="^70)
println("Example complete!")
println("="^70)

