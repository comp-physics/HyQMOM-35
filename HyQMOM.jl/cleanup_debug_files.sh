#!/bin/bash
# Cleanup temporary debugging files before final PR

echo "Cleaning up temporary test files..."

# Create archive directory
mkdir -p validation_archive

# Move one-time debugging tests to archive
mv -v test_rotational_invariance_final.jl validation_archive/ 2>/dev/null
mv -v test_unsplit_realizability.jl validation_archive/ 2>/dev/null
mv -v test_unsplit_physics.jl validation_archive/ 2>/dev/null
mv -v compare_split_unsplit_rotation.jl validation_archive/ 2>/dev/null

echo ""
echo "Cleanup complete!"
echo ""
echo "Files to KEEP (core validation):"
echo "  - test_moment_rotation_unit.jl"
echo "  - test_physical_rotation.jl"
echo "  - test_resolution_convergence.jl"
echo "  - test_split_vs_unsplit_crossing.jl"
echo "  - test_split_vs_unsplit_gridaligned.jl"
echo ""
echo "Files moved to validation_archive/:"
ls -1 validation_archive/ 2>/dev/null

