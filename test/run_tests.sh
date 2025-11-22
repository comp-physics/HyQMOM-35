#!/bin/bash
# Main Test Runner for HyQMOM.jl
#
# This is the consolidated master test runner that combines functionality from:
# - run_all_tests.sh (master runner)
# - ci_test_local.sh (local CI simulation)
#
# Usage:
#   ./test/run_tests.sh              # Run all tests
#   ./test/run_tests.sh unit         # Run only unit tests
#   ./test/run_tests.sh integration  # Run only integration tests
#   ./test/run_tests.sh mpi          # Run only MPI tests
#   ./test/run_tests.sh --quick      # Run quick tests only
#   ./test/run_tests.sh --help       # Show help

set -e

# Colors for output
if [ -t 1 ]; then
    GREEN='\033[0;32m'
    RED='\033[0;31m'
    YELLOW='\033[1;33m'
    BLUE='\033[0;34m'
    NC='\033[0m' # No Color
else
    GREEN=''
    RED=''
    YELLOW=''
    BLUE=''
    NC=''
fi

# Go to HyQMOM.jl root
cd "$(dirname "$0")/.."

# Parse arguments
TEST_TYPE="all"
QUICK_MODE=false

show_help() {
    cat << EOF
HyQMOM.jl Test Runner

Usage:
    ./test/run_tests.sh [OPTIONS] [TEST_TYPE]

Test Types:
    all          Run all tests (default)
    unit         Run only unit tests
    integration  Run only integration tests (Julia vs MATLAB)
    mpi          Run only MPI consistency tests

Options:
    --quick      Use quick test configurations
    --help       Show this help message

Examples:
    ./test/run_tests.sh                    # Run all tests
    ./test/run_tests.sh unit               # Run unit tests only
    ./test/run_tests.sh integration        # Run integration tests
    ./test/run_tests.sh mpi --quick        # Run quick MPI tests
EOF
    exit 0
}

for arg in "$@"; do
    case $arg in
        --help|-h)
            show_help
            ;;
        --quick)
            QUICK_MODE=true
            ;;
        unit|integration|mpi|all)
            TEST_TYPE=$arg
            ;;
        *)
            echo -e "${RED}Unknown argument: $arg${NC}"
            echo "Run with --help for usage information"
            exit 1
            ;;
    esac
done

echo -e "${BLUE}========================================================================"
echo "HYQMOM.JL - TEST SUITE"
echo "========================================================================"
echo -e "${NC}"
echo "Test type: $TEST_TYPE"
if [ "$QUICK_MODE" = true ]; then
    echo "Mode: Quick (reduced test configurations)"
fi
echo ""

# Check prerequisites
echo "Checking prerequisites..."

if ! command -v julia &> /dev/null; then
    echo -e "${RED}[X] Julia not found. Please install Julia 1.9 or later.${NC}"
    exit 1
fi
echo -e "${GREEN}[OK]${NC} Julia found: $(julia --version)"

if ! command -v mpiexec &> /dev/null && [ "$TEST_TYPE" = "mpi" -o "$TEST_TYPE" = "all" ]; then
    echo -e "${YELLOW}[WARNING]${NC} MPI not found. MPI tests will be skipped."
    SKIP_MPI=true
else
    SKIP_MPI=false
    if [ "$TEST_TYPE" = "mpi" -o "$TEST_TYPE" = "all" ]; then
        echo -e "${GREEN}[OK]${NC} MPI found: $(mpiexec --version | head -n 1)"
    fi
fi

# Check for MATLAB golden file (needed for integration tests)
GOLDEN_FILE="../goldenfiles/goldenfile_mpi_1ranks_Np20_tmax100.mat"
if [ "$TEST_TYPE" = "integration" -o "$TEST_TYPE" = "all" ]; then
    if [ ! -f "$GOLDEN_FILE" ]; then
        echo -e "${YELLOW}[WARNING]${NC} MATLAB golden file not found: $GOLDEN_FILE"
        echo "  Integration tests will be skipped."
        echo "  Run create_goldenfiles('ci') in MATLAB to generate it."
        SKIP_INTEGRATION=true
    else
        echo -e "${GREEN}[OK]${NC} MATLAB golden file found"
        SKIP_INTEGRATION=false
    fi
else
    SKIP_INTEGRATION=false
fi

echo ""

# Track failures
FAILED=0

# Function to print test result
print_result() {
    local exit_code=$1
    local test_name=$2
    
    if [ $exit_code -eq 0 ]; then
        echo -e "${GREEN}[OK] $test_name PASSED${NC}"
        return 0
    else
        echo -e "${RED}[X] $test_name FAILED (exit code: $exit_code)${NC}"
        return 1
    fi
}

# Step 1: Unit Tests
if [ "$TEST_TYPE" = "all" -o "$TEST_TYPE" = "unit" ]; then
    echo "========================================================================"
    echo "Step 1: Unit Tests"
    echo "========================================================================"
    
    # Skip integration tests in this phase by setting environment variable
    if TEST_INTEGRATION=false julia --project=. --color=yes -e 'using Pkg; Pkg.test()'; then
        print_result 0 "Unit tests"
    else
        print_result $? "Unit tests"
        FAILED=1
    fi
    echo ""
fi

# Step 2: Integration Tests (Julia vs MATLAB)
if [ "$TEST_TYPE" = "all" -o "$TEST_TYPE" = "integration" ]; then
    if [ "$SKIP_INTEGRATION" = true ]; then
        echo "========================================================================"
        echo "Step 2: Integration Tests - SKIPPED (no golden file)"
        echo "========================================================================"
        echo ""
    else
        echo "========================================================================"
        echo "Step 2: Integration Tests (Julia vs MATLAB)"
        echo "========================================================================"
        
        if julia --project=. test/test_integration.jl; then
            print_result 0 "Integration tests"
        else
            print_result $? "Integration tests"
            FAILED=1
        fi
        echo ""
    fi
fi

# Step 3: MPI Tests
if [ "$TEST_TYPE" = "all" -o "$TEST_TYPE" = "mpi" ]; then
    if [ "$SKIP_MPI" = true ]; then
        echo "========================================================================"
        echo "Step 3: MPI Tests - SKIPPED (MPI not available)"
        echo "========================================================================"
        echo ""
    else
        echo "========================================================================"
        echo "Step 3: MPI Consistency Tests"
        echo "========================================================================"
        
        if [ "$QUICK_MODE" = true ]; then
            # Quick mode: just test 1 vs 2 ranks
            echo "Quick mode: Testing 1 vs 2 ranks only"
            echo ""
            
            echo "Generating reference (1 rank)..."
            if julia --project=. test/test_mpi.jl; then
                echo -e "${GREEN}[OK] Reference generated${NC}"
            else
                echo -e "${RED}[X] Reference generation failed${NC}"
                FAILED=1
            fi
            
            echo ""
            echo "Testing with 2 ranks..."
            if mpiexec -n 2 julia --project=. test/test_mpi.jl; then
                print_result 0 "MPI tests (2 ranks)"
            else
                print_result $? "MPI tests (2 ranks)"
                FAILED=1
            fi
        else
            # Full mode: test multiple rank configurations
            if ./test/run_mpi_tests.sh; then
                print_result 0 "MPI tests"
            else
                print_result $? "MPI tests"
                FAILED=1
            fi
        fi
        echo ""
    fi
fi

# Summary
echo "========================================================================"
echo "TEST SUMMARY"
echo "========================================================================"

if [ $FAILED -eq 0 ]; then
    echo -e "${GREEN}"
    cat << EOF
+======================================================================+
|                    ALL TESTS PASSED! [OK][OK][OK]                            |
|                                                                      |
|  Your code is working correctly. Safe to commit and push.           |
+======================================================================+
EOF
    echo -e "${NC}"
    exit 0
else
    echo -e "${RED}"
    cat << EOF
+======================================================================+
|                    SOME TESTS FAILED [X][X][X]                            |
|                                                                      |
|  Please fix the failing tests before committing.                    |
+======================================================================+
EOF
    echo -e "${NC}"
    exit 1
fi

