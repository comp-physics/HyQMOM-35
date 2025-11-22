#!/usr/bin/env bash
# Restore Project.toml from git (if available) or add packages back
#
# Usage:
#   bash scripts/restore_project_toml.sh

set -euo pipefail

echo ""
echo "======================================================================"
echo "Restoring Project.toml"
echo "======================================================================"
echo ""

# Check if we're in a git repository
if git rev-parse --git-dir > /dev/null 2>&1; then
    echo "Git repository detected. Checking Project.toml status..."
    
    # Check if Project.toml has uncommitted changes
    if git diff --quiet HEAD -- Project.toml 2>/dev/null; then
        echo "✓ Project.toml is unchanged from git"
        echo "  No restoration needed."
    else
        echo "⚠ Project.toml has been modified"
        echo ""
        echo "Options:"
        echo "  1. Restore from git (discard local changes)"
        echo "  2. Keep current and add missing packages"
        echo "  3. Cancel"
        echo ""
        read -p "Choose option [1/2/3]: " choice
        
        case $choice in
            1)
                echo ""
                echo "Restoring Project.toml from git..."
                git checkout HEAD -- Project.toml
                echo "✓ Project.toml restored from git"
                echo ""
                echo "Next steps:"
                echo "  cd HyQMOM.jl"
                echo "  julia --project=. -e 'using Pkg; Pkg.instantiate()'"
                ;;
            2)
                echo ""
                echo "Adding missing packages..."
                julia --project=. scripts/restore_desktop.jl
                ;;
            3)
                echo "Cancelled."
                exit 0
                ;;
            *)
                echo "Invalid choice. Cancelled."
                exit 1
                ;;
        esac
    fi
else
    echo "Not a git repository. Adding packages via Pkg.add()..."
    julia --project=. scripts/restore_desktop.jl
fi

echo ""
echo "======================================================================"
echo "Done!"
echo "======================================================================"

