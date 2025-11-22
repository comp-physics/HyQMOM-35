# HyQMOM.jl

3D Hyperbolic Quadrature Method of Moments (HyQMOM) solver for the Boltzmann equation with BGK collision operator, featuring MPI parallelization and interactive visualization.

## Overview

HyQMOM.jl is a high-performance computational fluid dynamics solver that implements the Hyperbolic Quadrature Method of Moments (HyQMOM) for solving multidimensional kinetic equations. The solver extends the one-dimensional HyQMOM closure to handle 2D and 3D velocity distribution functions using up to 35 fourth-order velocity moments.

The method provides closures for higher-order moments without requiring quadrature points or velocity distribution function reconstruction. Instead, it uses orthogonal polynomials constructed directly from known moments to define the unclosed moments needed for the hyperbolic flux terms. This approach ensures global hyperbolicity and moment realizability across arbitrary Knudsen and Mach number regimes.

Key features include:

- **3D moment-based kinetic solver** for Boltzmann-BGK equation
- **MPI parallelization** with domain decomposition
- **Interactive 3D visualization** with GLMakie (time-series animation)
- **Flexible parameter control** via command-line arguments
- **Serial or parallel** execution (automatic detection)
- **Snapshot collection** for time evolution analysis

## Quick Navigation

```@contents
Pages = [
  "quickstart.md",
  "user_guide.md",
  "mpi.md",
  "tutorials/interactive_visualization.md",
  "tutorials/run_3d_jets_timeseries.md",
  "api.md",
  "dev_guide.md",
]
Depth = 2
```

## Getting Started

The fastest way to get started is with our Quickstart Guide, which will have you running your first simulation in minutes.

For a comprehensive overview of all features and capabilities, see the User Guide.

## Key Features

### High-Performance Computing
- MPI domain decomposition for scalable parallel execution
- Optimized numerical algorithms for kinetic equations
- Support for both serial and parallel workflows

### Interactive Visualization
- Real-time 3D isosurface visualization with GLMakie
- Time-series animation with interactive controls
- Multiple physical quantities (density, velocity, pressure, temperature)
- Works seamlessly with both serial and MPI parallel simulations

### Flexible Configuration
- Command-line parameter control for all simulation aspects
- Predefined initial condition configurations
- Easy customization for new physical scenarios

## License

MIT license.

See the [license file](https://github.com/comp-physics/HyQMOM.jl/blob/master/license.md) for the full text.

## Copyright

Copyright 2025 Spencer Bryngelson and Rodney Fox
