@echo off
REM HyQMOM.jl Documentation Build Script (Windows)
REM This script builds the documentation locally using Documenter.jl

echo 🔧 Building HyQMOM.jl Documentation...
echo ======================================

REM Check if we're in the right directory
if not exist "HyQMOM.jl\Project.toml" (
    echo ❌ Error: Please run this script from the repository root directory
    echo    Expected to find HyQMOM.jl\Project.toml
    exit /b 1
)

REM Navigate to HyQMOM.jl directory
cd HyQMOM.jl

echo 📦 Setting up documentation environment...

REM Set environment variable to skip heavy plotting dependencies
set HYQMOM_SKIP_PLOTTING=true

REM Install/update documentation dependencies
echo    - Installing documentation dependencies...
julia --project=docs -e "using Pkg; Pkg.develop(PackageSpec(path=\".\")); Pkg.instantiate()"

echo 📚 Building documentation...

REM Build the documentation
julia --project=docs docs/make.jl

if %ERRORLEVEL% equ 0 (
    echo.
    echo ✅ Documentation built successfully!
    echo.
    echo 📖 View the documentation:
    echo    Local file: %CD%\docs\build\html\index.html
    echo.
    echo 🌐 To serve locally ^(requires Python^):
    echo    cd docs\build\html ^&^& python -m http.server 8000
    echo    Then open: http://localhost:8000
    echo.
    echo 🔍 Documentation files:
    dir /b docs\build\html\*.html
) else (
    echo.
    echo ❌ Documentation build failed!
    echo    Check the error messages above for details.
    exit /b 1
)
