#!/bin/bash
# Test MPI consistency against golden files

set -e

cd "$(dirname "$0")/.."

echo "========================================================================"
echo "MPI GOLDEN FILE CONSISTENCY TESTS"
echo "========================================================================"
echo

# Check if golden files exist
GOLDEN_DIR="test/goldenfiles"
if [ ! -f "${GOLDEN_DIR}/mpi_1rank_small.bin" ]; then
    echo "ERROR: Golden files not found!"
    echo "Please run: ./test/create_mpi_goldenfiles.sh"
    exit 1
fi

# Test configurations
CONFIGS=("small" "medium")
RANKS=(2 4 8)

echo "Testing configurations: ${CONFIGS[@]}"
echo "Testing with ranks: ${RANKS[@]}"
echo

FAILED=0

for config in "${CONFIGS[@]}"; do
    echo "========================================================================"
    echo "Configuration: ${config}"
    echo "========================================================================"
    echo
    
    for nranks in "${RANKS[@]}"; do
        echo "------------------------------------------------------------------------"
        echo "Testing ${nranks} ranks vs 1 rank..."
        echo "------------------------------------------------------------------------"
        
        if mpiexec -n ${nranks} julia --project=. test/test_mpi_goldenfiles.jl ${config}; then
            echo "✓ ${nranks} ranks PASSED"
        else
            echo "✗ ${nranks} ranks FAILED"
            FAILED=1
        fi
        echo
    done
done

echo "========================================================================"
if [ $FAILED -eq 0 ]; then
    echo "✓ ALL MPI GOLDEN FILE TESTS PASSED"
else
    echo "✗ SOME TESTS FAILED"
fi
echo "========================================================================"

exit $FAILED

