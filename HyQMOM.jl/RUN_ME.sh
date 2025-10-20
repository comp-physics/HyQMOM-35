#!/bin/bash
# Quick start script for running simulation with standardized moment visualization

echo "Running simulation with standardized moments..."
echo "This will:"
echo "  1. Run crossing jets simulation"
echo "  2. Launch time-series viewer (press Enter when done)"
echo "  3. Save results to disk"
echo "  4. Launch scatterplot viewer for standardized moments"
echo ""

mpiexec -n 10 julia --project=. examples/run_3d_custom_jets.jl \
  --config crossing \
  --Ma 1.0 \
  --tmax 0.1 \
  --Nx 30 --Ny 30 --Nz 30 \
  --snapshot-interval 1 \
  --save-standardized-moments true

echo ""
echo "Done! Check for snapshots_*.jld2 file in current directory"

