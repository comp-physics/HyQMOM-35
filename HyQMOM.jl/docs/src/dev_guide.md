# Developer Guide

This guide covers the internal architecture of HyQMOM.jl and provides guidance for contributors and developers who want to extend or modify the codebase.

## Project Structure

### Directory Layout

```
HyQMOM.jl/
├── src/                          # Core source code
│   ├── HyQMOM.jl                # Main module file
│   ├── simulation_runner.jl     # Simulation orchestration
│   ├── initial_conditions.jl    # Initial condition setups
│   ├── autogen/                 # Auto-generated functions (MATLAB conversion)
│   │   ├── delta2star3D.jl
│   │   ├── jacobian6.jl
│   │   └── ...
│   ├── moments/                 # Moment-based operations
│   │   ├── InitializeM4_35.jl
│   │   ├── M2CS4_35.jl
│   │   ├── Moments5_3D.jl
│   │   ├── hyqmom_3D.jl
│   │   └── enforce_univariate.jl
│   ├── numerics/                # Numerical methods
│   │   ├── closure_and_eigenvalues.jl
│   │   ├── flux_HLL.jl
│   │   ├── collision35.jl
│   │   ├── apply_flux_update.jl
│   │   ├── apply_flux_update_3d.jl
│   │   └── ...
│   ├── realizability/           # Realizability constraints
│   │   ├── realizability.jl
│   │   ├── realizable_2D.jl
│   │   ├── realizable_3D.jl
│   │   └── ...
│   ├── mpi/                     # MPI parallelization
│   │   ├── setup_mpi_cartesian_3d.jl
│   │   ├── halo_exchange_3d.jl
│   │   ├── compute_halo_fluxes_and_wavespeeds_3d.jl
│   │   └── ...
│   ├── utils/                   # Utility functions
│   │   ├── moment_idx.jl
│   │   ├── compute_standardized_field.jl
│   │   ├── diagnostics.jl
│   │   └── ...
│   └── visualization/           # Visualization components
│       └── interactive_3d_timeseries_streaming.jl
├── examples/                    # Example scripts and configurations
│   ├── run_3d_jets_timeseries.jl
│   ├── run_3d_custom_jets.jl
│   └── parse_params.jl
├── test/                        # Test suite
├── docs/                        # Documentation source
├── visualize_jld2.jl           # Standalone visualization script
└── Project.toml                 # Package dependencies
```

### Core Components

**Main Simulation Pipeline:**
1. **Parameter parsing** (`examples/parse_params.jl`)
2. **Initial conditions** (`src/initial_conditions.jl`)
3. **Domain decomposition** (`src/mpi/setup_mpi_cartesian_3d.jl`)
4. **Time stepping loop** (`src/simulation_runner.jl`)
5. **Visualization/output** (`src/visualization/interactive_3d_timeseries_streaming.jl`)

**Numerical Core:**
- **Moment operations**: Converting between raw, central, and standardized moment representations
- **Realizability**: Ensuring physical validity of moment sets using matrix determinant constraints
- **HyQMOM closures**: Computing 5th-order moment closures from orthogonal polynomials
- **Flux computation**: Hyperbolic flux calculations with moment closures
- **Time integration**: Explicit time stepping schemes with realizability correction

## Architecture Principles

### Modular Design

HyQMOM.jl follows a modular architecture where each component has a specific responsibility:

- **Separation of concerns**: Numerics, MPI, visualization are separate modules
- **Composability**: Components can be mixed and matched for different use cases
- **Extensibility**: New initial conditions, numerical methods, or visualization modes can be added easily

### Performance Considerations

**Memory Layout:**
- 4D arrays: `[x, y, z, moment_component]`
- Contiguous memory access patterns for cache efficiency
- Minimal memory allocations in time-stepping loops

**Computational Efficiency:**
- In-place operations where possible
- Vectorized operations using Julia's broadcasting
- Specialized functions for hot code paths

### MPI Integration

**Design Philosophy:**
- MPI is integrated throughout, not bolted on
- All core functions are MPI-aware
- Graceful fallback to serial execution

**Communication Patterns:**
- Halo exchange for ghost cell updates
- Collective operations for global reductions
- Point-to-point communication for neighbor updates

## Adding New Features

### New Initial Conditions

To add a new initial condition configuration:

1. **Define the configuration** in `examples/parse_params.jl`:
   ```julia
   elseif config == "my_new_config"
       # Define jet parameters, domain setup, etc.
   ```

2. **Implement the initial condition** in `src/initial_conditions.jl`:
   ```julia
   function setup_my_new_initial_condition(params, grid)
       # Initialize moment fields
       return M_initial
   end
   ```

3. **Add documentation** and examples showing how to use it

### New Numerical Methods

To add a new flux method or time integrator:

1. **Create the implementation** in `src/numerics/`:
   ```julia
   function my_new_flux_method(M_left, M_right, params)
       # Compute numerical flux
       return flux
   end
   ```

2. **Integrate with the main loop** in `src/simulation_runner.jl`

3. **Add tests** in `test/` to verify correctness

4. **Update documentation** with usage examples

### New Visualization Features

To add new visualization capabilities:

1. **Create visualization function** in `src/visualization/`:
   ```julia
   function my_new_visualization(data, grid, params)
       # Implement visualization logic
   end
   ```

2. **Handle MPI considerations** (rank 0 only, data gathering, etc.)

3. **Add to the main visualization pipeline**

4. **Document the new feature** with examples

## Code Style and Conventions

### Naming Conventions

- **Functions**: `snake_case` (e.g., `compute_flux`, `update_moments`)
- **Types**: `PascalCase` (e.g., `GridParameters`, `SimulationState`)
- **Constants**: `UPPER_CASE` (e.g., `DEFAULT_CFL`, `MAX_ITERATIONS`)
- **Variables**: `snake_case` (e.g., `time_step`, `grid_size`)

### Documentation Standards

**Function Documentation:**
```julia
"""
    compute_flux(M_left, M_right, params)

Compute numerical flux between left and right states using HLL method.

# Arguments
- `M_left`: Left state moment vector
- `M_right`: Right state moment vector  
- `params`: Simulation parameters

# Returns
- `flux`: Numerical flux vector

# Examples
```julia
flux = compute_flux(M_L, M_R, params)
```

See also: [`flux_HLL`](@ref), [`Flux_closure35_and_realizable_3D`](@ref)
"""
function compute_flux(M_left, M_right, params)
    # Implementation
end
```

**Use DocStringExtensions for consistency:**
```julia
using DocStringExtensions

"""
$(SIGNATURES)

Brief description of the function.

$(FIELDS)  # For types
"""
```

### Testing Guidelines

**Test Organization:**
- Unit tests for individual functions
- Integration tests for complete workflows
- MPI tests for parallel functionality
- Visualization tests (when possible without display)

**Test Structure:**
```julia
@testset "My Feature Tests" begin
    @testset "Unit Tests" begin
        # Test individual functions
    end
    
    @testset "Integration Tests" begin
        # Test complete workflows
    end
    
    @testset "MPI Tests" begin
        # Test parallel functionality
    end
end
```

## Performance Optimization

### Profiling

**Built-in Julia profiler:**
```julia
using Profile
@profile run_simulation(params)
Profile.print()
```

**Memory profiling:**
```julia
using ProfileView
@profview run_simulation(params)
```

### Common Optimization Patterns

**Avoid allocations in hot loops:**
```julia
# Bad: allocates new arrays
function bad_update!(M, flux)
    M .= M + dt * flux  # Creates temporary array
end

# Good: in-place operations
function good_update!(M, flux, dt)
    @. M += dt * flux  # In-place update
end
```

**Use views for array slicing:**
```julia
# Bad: copies data
subarray = M[i:j, k:l, :, :]

# Good: creates view
subarray = @view M[i:j, k:l, :, :]
```

**Preallocate temporary arrays:**
```julia
struct SimulationCache
    flux_x::Array{Float64, 4}
    flux_y::Array{Float64, 4}
    temp_M::Array{Float64, 4}
end

function create_cache(grid_size)
    return SimulationCache(
        zeros(grid_size...),
        zeros(grid_size...),
        zeros(grid_size...)
    )
end
```

## MPI Development

### Adding MPI-Aware Functions

**Template for MPI-aware functions:**
```julia
function my_mpi_function(data, comm=MPI.COMM_WORLD)
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)
    
    if size == 1
        # Serial implementation
        return serial_version(data)
    else
        # Parallel implementation
        return parallel_version(data, rank, size, comm)
    end
end
```

**Halo exchange pattern:**
```julia
function exchange_halos!(M, neighbors, comm)
    # Send to neighbors, receive from neighbors
    for (direction, neighbor_rank) in neighbors
        if neighbor_rank != MPI.MPI_PROC_NULL
            # Extract halo data
            halo_data = extract_halo(M, direction)
            
            # Non-blocking send/receive
            req_send = MPI.Isend(halo_data, neighbor_rank, tag, comm)
            req_recv = MPI.Irecv!(recv_buffer, neighbor_rank, tag, comm)
            
            # Wait for completion
            MPI.Wait!(req_send)
            MPI.Wait!(req_recv)
            
            # Update ghost cells
            update_ghost_cells!(M, recv_buffer, direction)
        end
    end
end
```

## Visualization Development

### GLMakie Integration

**Basic visualization function structure:**
```julia
function my_visualization(data, grid, params)
    # Skip if plotting disabled
    if get(ENV, "HYQMOM_SKIP_PLOTTING", "false") == "true"
        @info "Plotting disabled, skipping visualization"
        return
    end
    
    # Only on rank 0 for MPI
    if MPI.Comm_rank(MPI.COMM_WORLD) != 0
        return
    end
    
    # GLMakie visualization code
    using GLMakie
    
    fig = Figure()
    # ... visualization implementation
    
    display(fig)
end
```

**Interactive controls pattern:**
```julia
function create_interactive_controls(fig, data)
    # Time slider
    time_slider = Slider(fig[2, 1], range=1:length(data), startvalue=1)
    
    # Quantity selection
    quantity_menu = Menu(fig[2, 2], options=["density", "u_velocity", "v_velocity"])
    
    # Connect controls to visualization updates
    on(time_slider.value) do t
        update_visualization!(t)
    end
    
    on(quantity_menu.selection) do qty
        update_quantity!(qty)
    end
end
```

## Contributing Guidelines

### Development Workflow

1. **Fork the repository** on GitHub
2. **Create a feature branch**: `git checkout -b feature/my-new-feature`
3. **Make changes** following the style guidelines
4. **Add tests** for new functionality
5. **Update documentation** as needed
6. **Run the test suite**: `julia --project=. -e 'using Pkg; Pkg.test()'`
7. **Submit a pull request** with a clear description

### Pull Request Guidelines

**Good pull request characteristics:**
- Clear, descriptive title
- Detailed description of changes
- Tests for new functionality
- Documentation updates
- No breaking changes (or clearly marked)
- Passes all existing tests

**Code review process:**
- Maintainers will review for correctness, style, and performance
- Address feedback promptly
- Squash commits before merging (if requested)

### Issue Reporting

**Good bug reports include:**
- Minimal reproducible example
- System information (OS, Julia version, MPI implementation)
- Error messages and stack traces
- Expected vs. actual behavior

**Feature requests should include:**
- Clear use case description
- Proposed API (if applicable)
- Willingness to contribute implementation

## Release Process

### Version Management

HyQMOM.jl follows semantic versioning (SemVer):
- **Major version** (X.0.0): Breaking changes
- **Minor version** (0.X.0): New features, backward compatible
- **Patch version** (0.0.X): Bug fixes, backward compatible

### Testing Requirements

Before release:
- All tests pass on supported Julia versions (1.9, 1.10, 1.11)
- MPI tests pass with different MPI implementations
- Documentation builds successfully
- Examples run without errors

### Documentation Updates

- Update version number in `Project.toml`
- Regenerate documentation with new version
- Tag release on GitHub
