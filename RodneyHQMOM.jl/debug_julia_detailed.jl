using RodneyHQMOM
using MPI
using Printf

# Initialize MPI
MPI.Init()

# Run simulation with detailed output
results = run_simulation(
    Np = 20,
    tmax = 0.1,
    num_workers = 1,
    verbose = false,  # We'll add our own output
    Kn = 1.0,
    Ma = 0.0,
    flag2D = 0,
    CFL = 0.5
)

# Print results
println("\nSimulation completed $(results[:nn]) steps")
println("Final M(1,1,1:5) = [$(results[:M][1,1,1]), $(results[:M][1,1,2]), $(results[:M][1,1,3]), $(results[:M][1,1,4]), $(results[:M][1,1,5])]")

MPI.Finalize()
