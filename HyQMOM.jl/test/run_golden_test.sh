#!/bin/bash
# Run golden file test comparing Julia vs MATLAB
#
# Usage:
#   ./run_golden_test.sh          # Run with Julia directly
#   ./run_golden_test.sh --mpi    # Run with MPI (1 rank)

set -e

cd "$(dirname "$0")/.."  # Go to HyQMOM.jl root

echo "========================================================================"
echo "Golden File Test: Julia vs MATLAB"
echo "========================================================================"
echo ""

if [[ "$1" == "--mpi" ]]; then
    echo "Running with MPI (1 rank)..."
    mpiexec -n 1 julia --project=. test/test_matlab_golden.jl
else
    echo "Running with Julia (no MPI)..."
    julia --project=. test/test_matlab_golden.jl
fi

exit_code=$?

echo ""
if [ $exit_code -eq 0 ]; then
    echo "========================================================================"
    echo "OK Golden file test PASSED"
    echo "========================================================================"
else
    echo "========================================================================"
    echo "FAIL Golden file test FAILED (exit code: $exit_code)"
    echo "========================================================================"
fi

exit $exit_code

