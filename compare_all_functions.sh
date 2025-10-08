#!/bin/bash
# Systematic comparison of all MATLAB and Julia functions

echo "=== SYSTEMATIC FUNCTION COMPARISON ==="
echo ""

# List of critical functions to compare
functions=(
    "collision35"
    "pas_HLL"
    "flux_HLL"
    "apply_flux_update"
    "closure_and_eigenvalues"
    "Flux_closure35_and_realizable_3D"
    "eigenvalues6_hyperbolic_3D"
    "compute_jacobian_eigenvalues"
    "InitializeM4_35"
    "M2CS4_35"
    "hyqmom_3D"
    "realizable_3D"
    "realizable_2D"
    "delta2star3D_permutation"
    "simulation_runner"
)

for func in "${functions[@]}"; do
    echo "----------------------------------------"
    echo "Function: $func"
    echo "----------------------------------------"
    
    # Find MATLAB version
    matlab_file=$(find src -name "${func}.m" 2>/dev/null | head -1)
    
    # Find Julia version
    julia_file=$(find RodneyHQMOM.jl/src -name "${func}.jl" 2>/dev/null | head -1)
    
    if [ -n "$matlab_file" ]; then
        echo "MATLAB: $matlab_file ($(wc -l < "$matlab_file") lines)"
    else
        echo "MATLAB: NOT FOUND"
    fi
    
    if [ -n "$julia_file" ]; then
        echo "Julia:  $julia_file ($(wc -l < "$julia_file") lines)"
    else
        echo "Julia:  NOT FOUND"
    fi
    
    echo ""
done

echo "========================================"
echo "SUMMARY: Compare line counts above"
echo "========================================"
