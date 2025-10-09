#!/bin/bash
# Comprehensive MPI Testing Suite
#
# Tests MPI parallelization by comparing results from different rank counts.
# Validates that domain decomposition produces identical results.
#
# Usage:
#   ./run_mpi_tests.sh              # Test 1 vs 2 ranks
#   ./run_mpi_tests.sh --extended   # Test 1 vs 2 vs 4 ranks

set -e

cd "$(dirname "$0")/.."  # Go to HyQMOM.jl root

EXTENDED=false
if [[ "$1" == "--extended" ]]; then
    EXTENDED=true
fi

echo "========================================================================"
echo "MPI PARALLEL TESTING SUITE FOR HyQMOM.jl"
echo "========================================================================"
echo ""

# Clean up any previous test artifacts
echo "Cleaning up previous test results..."
rm -f test/mpi_reference_*.bin test/mpi_results_*.bin

# Test 1: Generate reference with 1 rank
echo ""
echo "========================================================================"
echo "TEST 1: Generate reference data (1 rank)"
echo "========================================================================"
julia --project=. test/test_mpi_integration.jl
EXIT_1=$?

if [ $EXIT_1 -ne 0 ]; then
    echo ""
    echo "FAIL: 1-rank test failed (exit code: $EXIT_1)"
    exit $EXIT_1
fi

# Test 2: Compare with 2 ranks
echo ""
echo "========================================================================"
echo "TEST 2: Verify MPI consistency (2 ranks)"
echo "========================================================================"
mpiexec -n 2 julia --project=. test/test_mpi_integration.jl
EXIT_2=$?

if [ $EXIT_2 -ne 0 ]; then
    echo ""
    echo "FAIL: 2-rank test failed (exit code: $EXIT_2)"
    exit $EXIT_2
fi

# Test 3: Extended test with 4 ranks (optional)
if [ "$EXTENDED" = true ]; then
    echo ""
    echo "========================================================================"
    echo "TEST 3: Extended MPI consistency test (4 ranks)"
    echo "========================================================================"
    mpiexec -n 4 julia --project=. test/test_mpi_integration.jl
    EXIT_4=$?
    
    if [ $EXIT_4 -ne 0 ]; then
        echo ""
        echo "FAIL: 4-rank test failed (exit code: $EXIT_4)"
        exit $EXIT_4
    fi
fi

# All tests passed
echo ""
echo "========================================================================"
echo "OK: ALL MPI TESTS PASSED"
echo "========================================================================"
echo ""
echo "Summary:"
echo "  ✓ 1-rank simulation runs successfully"
echo "  ✓ 2-rank simulation produces identical results"
if [ "$EXTENDED" = true ]; then
    echo "  ✓ 4-rank simulation produces identical results"
fi
echo ""
echo "MPI parallelization is working correctly!"
echo "========================================================================"

# Clean up test artifacts
echo ""
echo "Cleaning up test files..."
rm -f test/mpi_reference_*.bin test/mpi_results_*.bin

exit 0

