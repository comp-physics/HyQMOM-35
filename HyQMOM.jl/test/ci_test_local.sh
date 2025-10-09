#!/bin/bash
# Local CI Test Script
# Simulates what GitHub Actions will run
#
# Usage:
#   ./ci_test_local.sh           # Run all tests
#   ./ci_test_local.sh unit      # Run only unit tests
#   ./ci_test_local.sh golden    # Run only golden file test
#   ./ci_test_local.sh mpi       # Run only MPI test

set -e

cd "$(dirname "$0")/.."  # Go to HyQMOM.jl root

TEST_TYPE="${1:-all}"

echo "========================================================================"
echo "Local CI Test Simulation"
echo "========================================================================"
echo "Test type: $TEST_TYPE"
echo ""

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print status
print_status() {
    if [ $1 -eq 0 ]; then
        echo -e "${GREEN}OK $2 PASSED${NC}"
    else
        echo -e "${RED}FAIL $2 FAILED${NC}"
        return 1
    fi
}

# Check prerequisites
echo "Checking prerequisites..."
if ! command -v julia &> /dev/null; then
    echo -e "${RED}FAIL Julia not found. Please install Julia 1.9 or later.${NC}"
    exit 1
fi

if ! command -v mpiexec &> /dev/null && [ "$TEST_TYPE" = "mpi" -o "$TEST_TYPE" = "all" ]; then
    echo -e "${YELLOW}WARNING  MPI not found. MPI tests will be skipped.${NC}"
    SKIP_MPI=true
else
    SKIP_MPI=false
fi

# Check golden file
GOLDEN_FILE="../goldenfiles/goldenfile_mpi_1ranks_Np20_tmax100.mat"
if [ ! -f "$GOLDEN_FILE" ]; then
    echo -e "${RED}FAIL Golden file not found: $GOLDEN_FILE${NC}"
    echo "   Run create_goldenfiles('ci') in MATLAB to generate it."
    exit 1
fi
echo -e "${GREEN}OK${NC} Golden file found"

# Check Julia version
JULIA_VERSION=$(julia --version | grep -oE '[0-9]+\.[0-9]+')
echo -e "${GREEN}OK${NC} Julia version: $JULIA_VERSION"
echo ""

# Test 1: Install dependencies
if [ "$TEST_TYPE" = "all" -o "$TEST_TYPE" = "unit" ]; then
    echo "========================================================================"
    echo "Step 1: Installing Julia dependencies"
    echo "========================================================================"
    julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.status()'
    print_status $? "Dependency installation"
    echo ""
fi

# Test 2: Build MPI
if [ "$TEST_TYPE" = "all" -o "$TEST_TYPE" = "mpi" ] && [ "$SKIP_MPI" = false ]; then
    echo "========================================================================"
    echo "Step 2: Building MPI.jl"
    echo "========================================================================"
    julia --project=. -e 'using Pkg; Pkg.build("MPI"); using MPI; println("MPI OK")'
    print_status $? "MPI build"
    echo ""
fi

# Test 3: Unit tests
if [ "$TEST_TYPE" = "all" -o "$TEST_TYPE" = "unit" ]; then
    echo "========================================================================"
    echo "Step 3: Running Julia unit tests"
    echo "========================================================================"
    # Set environment variable to skip golden test in Pkg.test() since we run it separately
    TEST_GOLDEN_SIMULATION=false julia --project=. --color=yes -e 'using Pkg; Pkg.test()'
    UNIT_TEST_STATUS=$?
    print_status $UNIT_TEST_STATUS "Unit tests"
    echo ""
fi

# Test 4: Golden file test
if [ "$TEST_TYPE" = "all" -o "$TEST_TYPE" = "golden" ]; then
    echo "========================================================================"
    echo "Step 4: Running golden file test (Julia vs MATLAB)"
    echo "========================================================================"
    julia --project=. test/test_matlab_golden.jl
    GOLDEN_TEST_STATUS=$?
    print_status $GOLDEN_TEST_STATUS "Golden file test"
    echo ""
fi

# Test 5: MPI golden file test
if [ "$TEST_TYPE" = "all" -o "$TEST_TYPE" = "mpi" ] && [ "$SKIP_MPI" = false ]; then
    echo "========================================================================"
    echo "Step 5: Running MPI golden file test (1 rank)"
    echo "========================================================================"
    mpiexec -n 1 julia --project=. test/test_matlab_golden.jl
    MPI_TEST_STATUS=$?
    print_status $MPI_TEST_STATUS "MPI golden file test"
    echo ""
fi

# Summary
echo "========================================================================"
echo "LOCAL CI TEST SUMMARY"
echo "========================================================================"

EXIT_CODE=0

if [ "$TEST_TYPE" = "all" -o "$TEST_TYPE" = "unit" ]; then
    if [ ${UNIT_TEST_STATUS:-0} -eq 0 ]; then
        echo -e "${GREEN}OK Unit tests: PASSED${NC}"
    else
        echo -e "${RED}FAIL Unit tests: FAILED${NC}"
        EXIT_CODE=1
    fi
fi

if [ "$TEST_TYPE" = "all" -o "$TEST_TYPE" = "golden" ]; then
    if [ ${GOLDEN_TEST_STATUS:-0} -eq 0 ]; then
        echo -e "${GREEN}OK Golden file test: PASSED${NC}"
    else
        echo -e "${RED}FAIL Golden file test: FAILED${NC}"
        EXIT_CODE=1
    fi
fi

if [ "$TEST_TYPE" = "all" -o "$TEST_TYPE" = "mpi" ] && [ "$SKIP_MPI" = false ]; then
    if [ ${MPI_TEST_STATUS:-0} -eq 0 ]; then
        echo -e "${GREEN}OK MPI golden file test: PASSED${NC}"
    else
        echo -e "${RED}FAIL MPI golden file test: FAILED${NC}"
        EXIT_CODE=1
    fi
fi

echo "========================================================================"

if [ $EXIT_CODE -eq 0 ]; then
    echo -e "${GREEN}"
    echo "+==================================================================+"
    echo "|                 ALL TESTS PASSED! SUCCESS                             |"
    echo "|                                                                  |"
    echo "|  Your code is ready for CI. GitHub Actions will run these       |"
    echo "|  same tests automatically on push/PR.                           |"
    echo "+==================================================================+"
    echo -e "${NC}"
else
    echo -e "${RED}"
    echo "+==================================================================+"
    echo "|                   SOME TESTS FAILED                              |"
    echo "|                                                                  |"
    echo "|  Please fix the failing tests before pushing to GitHub.         |"
    echo "|  CI will run the same tests and fail if these aren't fixed.     |"
    echo "+==================================================================+"
    echo -e "${NC}"
fi

exit $EXIT_CODE

