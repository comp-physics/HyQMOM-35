#!/bin/bash
# MPI Testing Suite
#
# Tests MPI parallelization by comparing results from different rank counts.
# This script now uses the consolidated test/test_mpi.jl
#
# Usage:
#   ./run_mpi_tests.sh              # Test 1 vs 2 ranks
#   ./run_mpi_tests.sh --extended   # Test 1 vs 2 vs 4 ranks
#   ./run_mpi_tests.sh --golden     # Test against pre-generated golden files

set -e

cd "$(dirname "$0")/.."  # Go to HyQMOM.jl root

EXTENDED=false
GOLDEN_MODE=false

for arg in "$@"; do
    case $arg in
        --extended)
            EXTENDED=true
            ;;
        --golden)
            GOLDEN_MODE=true
            ;;
        *)
            echo "Unknown option: $arg"
            echo "Usage: ./run_mpi_tests.sh [--extended] [--golden]"
            exit 1
            ;;
    esac
done

echo "========================================================================"
echo "MPI PARALLEL TESTING SUITE FOR HyQMOM.jl"
echo "========================================================================"
echo ""

# Clean up any previous test artifacts
echo "Cleaning up previous test results..."
rm -f test/mpi_reference_*.bin test/mpi_results_*.bin

if [ "$GOLDEN_MODE" = true ]; then
    # Golden file mode - test against pre-generated files
    echo "Mode: Golden file testing"
    echo ""
    
    # Check if golden files exist
    if [ ! -f "test/goldenfiles/mpi_1rank_small.bin" ]; then
        echo "ERROR: Golden files not found!"
        echo "Generate them with: mpiexec -n 1 julia --project=. test/create_golden_files.jl"
        exit 1
    fi
    
    # Test small configuration
    echo "========================================================================"
    echo "TEST 1: Small configuration (2 ranks)"
    echo "========================================================================"
    mpiexec -n 2 julia --project=. test/test_mpi.jl --golden small
    EXIT_1=$?
    
    if [ $EXIT_1 -ne 0 ]; then
        echo "FAIL: Small config test failed (exit code: $EXIT_1)"
        exit $EXIT_1
    fi
    
    # Test medium configuration if extended
    if [ "$EXTENDED" = true ]; then
        echo ""
        echo "========================================================================"
        echo "TEST 2: Medium configuration (4 ranks)"
        echo "========================================================================"
        mpiexec -n 4 julia --project=. test/test_mpi.jl --golden medium
        EXIT_2=$?
        
        if [ $EXIT_2 -ne 0 ]; then
            echo "FAIL: Medium config test failed (exit code: $EXIT_2)"
            exit $EXIT_2
        fi
    fi
    
else
    # Dynamic mode - generate reference and compare
    echo "Mode: Dynamic comparison"
    echo ""
    
    # Test 1: Generate reference with 1 rank
    echo "========================================================================"
    echo "TEST 1: Generate reference data (1 rank)"
    echo "========================================================================"
    julia --project=. test/test_mpi.jl
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
    mpiexec -n 2 julia --project=. test/test_mpi.jl
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
        mpiexec -n 4 julia --project=. test/test_mpi.jl
        EXIT_4=$?
        
        if [ $EXIT_4 -ne 0 ]; then
            echo ""
            echo "FAIL: 4-rank test failed (exit code: $EXIT_4)"
            exit $EXIT_4
        fi
    fi
fi

# All tests passed
echo ""
echo "========================================================================"
echo "✓ ALL MPI TESTS PASSED"
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

