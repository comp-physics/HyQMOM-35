#!/bin/bash
# Create MPI Golden Files for Regression Testing
#
# This script generates golden reference files for MPI consistency testing.
# It consolidates functionality from create_mpi_goldenfiles.sh
#
# Usage:
#   ./test/create_golden_files.sh              # Generate all configurations
#   ./test/create_golden_files.sh small        # Generate only small config
#   ./test/create_golden_files.sh small medium # Generate specific configs

set -e

cd "$(dirname "$0")/.."  # Go to HyQMOM.jl root

echo "========================================================================"
echo "CREATING MPI GOLDEN FILES"
echo "========================================================================"
echo ""

# Check for MPI
if ! command -v mpiexec &> /dev/null; then
    echo "ERROR: mpiexec not found. Please install MPI."
    exit 1
fi

# Run the Julia script to generate golden files
echo "Generating golden files with 1 MPI rank..."
echo ""

if [ $# -eq 0 ]; then
    # No arguments - generate all
    mpiexec -n 1 julia --project=. test/create_golden_files.jl
else
    # Pass arguments to Julia script
    mpiexec -n 1 julia --project=. test/create_golden_files.jl "$@"
fi

exit_code=$?

if [ $exit_code -eq 0 ]; then
    echo ""
    echo "========================================================================"
    echo "✓ GOLDEN FILES CREATED SUCCESSFULLY"
    echo "========================================================================"
    echo ""
    echo "Golden files saved in: test/goldenfiles/"
    echo ""
    echo "To test MPI consistency:"
    echo "  ./test/run_mpi_tests.sh --golden"
    echo ""
else
    echo ""
    echo "========================================================================"
    echo "✗ GOLDEN FILE GENERATION FAILED"
    echo "========================================================================"
    exit $exit_code
fi

