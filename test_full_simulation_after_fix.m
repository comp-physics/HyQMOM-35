% Test full simulation after bug fix

clear all;
close all;

setup_paths;

fprintf('===== TESTING FULL SIMULATION AFTER FIX - MATLAB =====\n\n');

params = struct('Np', 20, 'Kn', 1.0, 'Ma', 0.0, 'tmax', 0.1, 'flag2D', 0, 'CFL', 0.5);

fprintf('Running simulation with tmax = 0.1...\n');
try
    result = simulation_runner(params);
    
    fprintf('\n=== SIMULATION COMPLETED ===\n');
    fprintf('Time steps: %d\n', result.time_steps);
    fprintf('Final time: %.6f\n', result.final_time);
    fprintf('Steps completed: %d\n', result.time_steps);
    
    % Check for NaN
    if any(isnan(result.M(:)))
        fprintf('❌ FAILED: NaN detected in final moments\n');
    else
        fprintf('✅ SUCCESS: No NaN in simulation\n');
    end
    
    % Check moment values
    M_final = result.M;
    max_moment = max(abs(M_final(:)));
    fprintf('Max moment magnitude: %.6e\n', max_moment);
    
    if max_moment > 1e6
        fprintf('❌ WARNING: Moments exploded\n');
    else
        fprintf('✅ Moments stayed bounded\n');
    end
    
catch e
    fprintf('❌ SIMULATION FAILED:\n');
    fprintf('%s\n', e.message);
end

fprintf('\n===== END =====\n');
