#!/bin/bash
# Run MPI consistency test: 1 rank vs 2 ranks
#
# This script runs the same simulation with 1 and 2 MPI ranks,
# then compares the results to verify MPI consistency.

set -e

cd "$(dirname "$0")"

# Clean up any previous results
echo "========================================================================"
echo "MPI CONSISTENCY TEST SUITE"
echo "========================================================================"
echo ""
echo "Cleaning up previous test results..."
rm -f mpi_results_1ranks_*.bin mpi_results_2ranks_*.bin

echo ""
echo "========================================================================"
echo "Step 1: Running with 1 MPI rank"
echo "========================================================================"
julia --project=.. test_mpi_consistency.jl

exit_code_1=$?

if [ $exit_code_1 -ne 0 ]; then
    echo ""
    echo "========================================================================"
    echo "FAIL: 1-rank simulation failed (exit code: $exit_code_1)"
    echo "========================================================================"
    exit $exit_code_1
fi

echo ""
echo "========================================================================"
echo "Step 2: Running with 2 MPI ranks"
echo "========================================================================"
mpiexec -n 2 julia --project=.. test_mpi_consistency.jl

exit_code_2=$?

echo ""
echo "========================================================================"
if [ $exit_code_2 -eq 0 ]; then
    echo "OK: MPI consistency test PASSED"
    echo "    Results from 1 and 2 ranks match!"
else
    echo "FAIL: MPI consistency test FAILED (exit code: $exit_code_2)"
    echo "    Results from 1 and 2 ranks differ!"
fi
echo "========================================================================"

# Clean up test artifacts
echo ""
echo "Cleaning up test files..."
rm -f mpi_results_1ranks_*.bin mpi_results_2ranks_*.bin

exit $exit_code_2

