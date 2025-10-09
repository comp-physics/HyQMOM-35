#!/bin/bash
# Create MPI golden files for regression testing

set -e

cd "$(dirname "$0")/.."

echo "========================================================================"
echo "CREATING MPI GOLDEN FILES"
echo "========================================================================"
echo

# Create golden files with 1 rank
echo "Generating 1-rank reference golden files..."
mpiexec -n 1 julia --project=. test/create_mpi_goldenfiles.jl

echo
echo "========================================================================"
echo "GOLDEN FILES CREATED"
echo "========================================================================"
echo
echo "To test MPI consistency, run:"
echo "  ./test/test_mpi_goldenfiles.sh"

