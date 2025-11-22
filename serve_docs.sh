#!/bin/bash

# HyQMOM.jl Documentation Server Script
# This script serves the built documentation locally

set -e

echo "🌐 Starting HyQMOM.jl Documentation Server..."
echo "============================================="

# Check if documentation exists
if [ ! -f "docs/build/html/index.html" ]; then
    echo "❌ Error: Documentation not found!"
    echo "   Please build the documentation first:"
    echo "   ./build_docs.sh"
    exit 1
fi

# Navigate to documentation directory
cd docs/build/html

echo "📖 Documentation server starting..."
echo "   URL: http://localhost:8000"
echo "   Press Ctrl+C to stop the server"
echo ""

# Try different Python commands
if command -v python3 &> /dev/null; then
    python3 -m http.server 8000
elif command -v python &> /dev/null; then
    python -m http.server 8000
else
    echo "❌ Error: Python not found!"
    echo "   Please install Python to serve the documentation"
    echo "   Alternatively, open the file directly:"
    echo "   file://$(pwd)/index.html"
    exit 1
fi
