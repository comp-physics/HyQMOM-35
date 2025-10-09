#!/bin/bash
# Master test runner - runs all tests in the proper order
#
# This is the main entry point for CI and local testing.
#
# Usage:
#   ./test/run_all_tests.sh

set -e

cd "$(dirname "$0")/.."  # Go to HyQMOM.jl root

echo "========================================================================"
echo "HYQMOM.JL - COMPLETE TEST SUITE"
echo "========================================================================"
echo

FAILED=0

# Step 1: Unit tests
echo "========================================================================"
echo "Step 1: Unit Tests"
echo "========================================================================"
if TEST_GOLDEN_SIMULATION=false CI=true julia --project=. --color=yes -e 'using Pkg; Pkg.test()'; then
    echo "✓ Unit tests PASSED"
else
    echo "✗ Unit tests FAILED"
    FAILED=1
fi
echo

# Step 2: MPI consistency tests (golden files)
echo "========================================================================"
echo "Step 2: MPI Golden File Tests (Regression)"
echo "========================================================================"
if ./test/test_mpi_goldenfiles.sh; then
    echo "✓ MPI golden file tests PASSED"
else
    echo "✗ MPI golden file tests FAILED"
    FAILED=1
fi
echo

# Step 3: MPI integration tests (dynamic)
echo "========================================================================"
echo "Step 3: MPI Integration Tests (Dynamic)"
echo "========================================================================"
if ./test/run_mpi_tests.sh; then
    echo "✓ MPI integration tests PASSED"
else
    echo "✗ MPI integration tests FAILED"
    FAILED=1
fi
echo

# Summary
echo "========================================================================"
if [ $FAILED -eq 0 ]; then
    echo "✓✓✓ ALL TESTS PASSED ✓✓✓"
    echo "========================================================================"
    exit 0
else
    echo "✗✗✗ SOME TESTS FAILED ✗✗✗"
    echo "========================================================================"
    exit 1
fi

