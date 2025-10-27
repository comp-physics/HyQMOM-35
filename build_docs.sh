#!/bin/bash

# HyQMOM.jl Documentation Build Script
# This script builds the documentation locally using Documenter.jl

set -e  # Exit on any error

echo "🔧 Building HyQMOM.jl Documentation..."
echo "======================================"

# Check if we're in the right directory
if [ ! -f "HyQMOM.jl/Project.toml" ]; then
    echo "❌ Error: Please run this script from the repository root directory"
    echo "   Expected to find HyQMOM.jl/Project.toml"
    exit 1
fi

# Navigate to HyQMOM.jl directory
cd HyQMOM.jl

echo "📦 Setting up documentation environment..."

# Set environment variable to skip heavy plotting dependencies
export HYQMOM_SKIP_PLOTTING=true

# Install/update documentation dependencies
echo "   - Installing documentation dependencies..."
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=".")); Pkg.instantiate()'

echo "📚 Building documentation..."

# Build the documentation
julia --project=docs docs/make.jl

# Check if build was successful
if [ $? -eq 0 ]; then
    echo ""
    echo "✅ Documentation built successfully!"
    echo ""
    echo "📖 View the documentation:"
    echo "   Local file: $(pwd)/docs/build/html/index.html"
    echo ""
    echo "🌐 To serve locally (requires Python):"
    echo "   cd HyQMOM.jl/docs/build/html && python3 -m http.server 4000"
    echo "   Or use the convenience script: ./serve_docs.sh"
    echo "   Then open: http://localhost:4000"
    echo ""
    echo "🔍 Documentation files:"
    ls -la docs/build/html/*.html | head -10
    if [ $(ls docs/build/html/*.html | wc -l) -gt 10 ]; then
        echo "   ... and $(( $(ls docs/build/html/*.html | wc -l) - 10 )) more files"
    fi
else
    echo ""
    echo "❌ Documentation build failed!"
    echo "   Check the error messages above for details."
    exit 1
fi
