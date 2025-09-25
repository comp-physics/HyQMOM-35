# 3D HyQMOM Riemann Solver

A MATLAB implementation of a 2D Riemann solver for 3D Hyperbolic Quadrature Method of Moments (HyQMOM) using the HLL (Harten-Lax-van Leer) scheme and explicit Euler time integration.

## Overview

This code solves the Boltzmann equation using moment methods in a 3D velocity space with 2D spatial discretization.
The implementation is restricted to 4th-order moments in 3D with 35 total moments, providing a high-fidelity simulation of kinetic gas dynamics with collision effects.

### Key Features

- **3D HyQMOM**: Hyperbolic Quadrature Method of Moments for kinetic theory
- **2D Spatial Domain**: Square computational domain with customizable grid resolution
- **35 Moment System**: Complete 4th-order moment closure in 3D velocity space
- **HLL Riemann Solver**: Robust flux computation for hyperbolic systems
- **BGK Collision Model**: Bhatnagar-Gross-Krook collision operator
- **Realizability Enforcement**: Ensures physical validity of moment solutions
- **Hyperbolicity Checks**: Maintains mathematical well-posedness

## Mathematical Background

The code solves the moment equations derived from the Boltzmann equation:

$$\frac{\partial \mathbf{M}}{\partial t} + \nabla \cdot \mathbf{F}(\mathbf{M}) = \mathbf{C}(\mathbf{M})$$

Where:
- $\mathbf{M}$ = 35-component moment vector
- $\mathbf{F}(\mathbf{M})$ = Flux tensor (computed via HyQMOM)
- $\mathbf{C}(\mathbf{M})$ = BGK collision operator

### Moment Vector Structure

The 35 moments are organized as $\mathbf{M} = [M_{000}, M_{100}, M_{200}, \ldots, M_{022}]$ representing:
- $M_{ijk} = \int u^i v^j w^k f(u,v,w) \, du \, dv \, dw$
- Density, velocities, temperatures, and higher-order moments

## Installation & Requirements

### Prerequisites
- MATLAB R2020b or later
- Parallel Computing Toolbox (optional, for performance)

## Usage

### Basic Simulation

```matlab
% Run with default parameters
results = main();

% Custom grid resolution and final time
results = main(50, 0.1);

% Disable plotting for faster execution
results = main(50, 0.1, false);

% Enable output file saving
results = main(50, 0.1, false, true);
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `Np` | 6 | Grid points per dimension (total: Np×Np) |
| `tmax` | 0.05 | Final simulation time |
| `enable_plots` | true | Generate visualization plots |
| `save_output` | false | Save results to .mat file |

### Example Workflows

**Quick Test Run:**
```matlab
results = main(6, 0.02, false);  % Small grid, short time, no plots
```

**Production Run:**
```matlab
results = main(100, 0.5, true, true);  % High resolution, full time, with plots and saving
```

**Batch Processing:**
```matlab
for Np = [10, 20, 50, 100]
    results = main(Np, 0.1, false, true);
    fprintf('Completed Np=%d\n', Np);
end
```

## Physical Problem Setup

The default configuration simulates a **crossing jets problem**:

- **Domain**: $[-0.5, 0.5] \times [-0.5, 0.5]$
- **Initial Conditions**: 
  - Low-density background ($\rho = 0.01$)
  - High-density crossing jets ($\rho = 1.0$)
  - Temperature $T = 1$ (dimensionless)
- **Boundary Conditions**: Periodic
- **Knudsen Number**: $\text{Kn} = 1$ (moderate rarefaction)
- **Mach Number**: $\text{Ma} = 0$ (subsonic)

### Physical Parameters

- **CFL Number**: 0.5 (stability constraint)
- **Collision Frequency**: $1/\text{Kn}$ (BGK model)
- **Gas Model**: Monatomic ideal gas

## Output Data

The `results` structure contains:

```matlab
results.parameters    % Simulation parameters (Np, tmax, Kn, Ma, etc.)
results.grid         % Spatial discretization (x, y, dx, dy)
results.moments      % Final moment fields (M, C, S, M5, C5, S5)
results.eigenvalues  % Characteristic speeds (lam6xa, lam6xb, etc.)
results.velocities   % Wave speed bounds (v5xmin, v5xmax, etc.)
results.filename     % Output file name (if saved)
```

## Visualization

The code generates several types of plots:

1. **Initial Conditions** (Figures 1-4)
   - Moment contours
   - Central moments
   - Standardized moments

2. **Time Evolution** (Figure 10)
   - Real-time moment evolution
   - Collision effects

3. **Final Results** (Figures 5-9)
   - Final moment distributions
   - Eigenvalue analysis
   - Realizability diagnostics

## Testing & Validation

### Golden File System

The code includes a robust testing framework:

```matlab
% Create reference solution
run_goldenfile_creation

% Validate against reference
runtests('tests/test_validate_goldenfile.m')
```

### Continuous Integration

GitHub Actions automatically validates all changes:
- Runs regression tests
- Compares against golden files
- Ensures numerical accuracy (tolerance: $10^{-10}$)

## Performance

### Computational Complexity
- **Spatial**: $O(N_p^2)$ grid points
- **Moments**: 35 coupled equations per grid point
- **Time Steps**: Adaptive (CFL-limited)

### Optimization Tips
- Use `enable_plots = false` for production runs
- Enable Parallel Computing Toolbox for large grids
- Consider `save_output = false` to reduce I/O

## Numerical Methods

### Spatial Discretization
- **Finite Volume**: Conservative formulation
- **HLL Riemann Solver**: Robust flux computation
- **Second-Order Accuracy**: Spatial reconstruction

### Time Integration
- **Explicit Euler**: First-order time stepping
- **Operator Splitting**: Advection + collision
- **CFL Stability**: Adaptive time step control

### Moment Closure
- **HyQMOM**: Hyperbolic quadrature method
- **Realizability**: Enforced at each time step
- **Eigenvalue Analysis**: Hyperbolicity verification