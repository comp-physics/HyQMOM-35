"""
Wrapper to run simulation with snapshot saving

This adds the snapshot_interval parameter to params and calls simulation_runner.
"""

"""
    run_simulation_with_snapshots(params; snapshot_interval=10)

Run simulation and collect snapshots at regular intervals.

This is a convenience wrapper that adds `snapshot_interval` to the params
and calls `simulation_runner`, which now has built-in snapshot support.

# Arguments
- `params`: Named tuple with simulation parameters
- `snapshot_interval`: Save snapshot every N steps (default: 10)

# Returns  
On rank 0:
- `snapshots`: Vector of (M, t, step) named tuples
- `grid`: Grid structure

On other ranks:
- `nothing, nothing`
"""
function run_simulation_with_snapshots(params; snapshot_interval=10)
    # Add snapshot_interval to params
    params_with_snapshots = merge(params, (snapshot_interval=snapshot_interval,))
    
    # Call simulation_runner with snapshot support
    return simulation_runner(params_with_snapshots)
end

"""
    simulation_runner_snapshots(params)

Alias for run_simulation_with_snapshots for backward compatibility.

Extracts snapshot_interval from params if present, otherwise uses default.
"""
function simulation_runner_snapshots(params)
    snapshot_interval = get(params, :snapshot_interval, 10)
    return run_simulation_with_snapshots(params; snapshot_interval=snapshot_interval)
end

