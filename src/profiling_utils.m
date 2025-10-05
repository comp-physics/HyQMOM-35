function varargout = profiling_utils(action, varargin)
% Utilities for MPI profiling management
% Handles initialization, data collection, and reporting
%
% Usage:
%   profiling_utils('start')                    - Initialize MPI profiler
%   profiling_utils('report', num_workers, Np)  - Collect and display results

switch action
    case 'start'
        start_profiling();
        
    case 'report'
        if nargin < 3
            error('profiling_utils(''report'') requires num_workers and Np arguments');
        end
        num_workers = varargin{1};
        Np = varargin{2};
        report_profiling(num_workers, Np);
        
    otherwise
        error('Unknown action: %s', action);
end

end

function start_profiling()
% Initialize MPI profiler
    fprintf('MPI profiling enabled. Starting mpiprofile...\n');
    mpiprofile on;
end

function report_profiling(num_workers, Np)
% Collect and display MPI profiling results
    fprintf('\nCollecting MPI profile data...\n');
    profile_stats = mpiprofile('info');
    
    % Save profile data to file
    profile_filename = sprintf('mpi_profile_%dranks_Np%d.mat', num_workers, Np);
    save(profile_filename, 'profile_stats');
    fprintf('Profile data saved to: %s\n', profile_filename);
    fprintf('To view later, use: load(''%s''); mpiprofile(''viewer'', profile_stats)\n', profile_filename);
    
    % Open profile viewer
    fprintf('Opening MPI profile viewer...\n');
    mpiprofile viewer;
    
    % Display summary statistics
    fprintf('\n=== MPI Profile Summary ===\n');
    for w = 1:length(profile_stats)
        if ~isempty(profile_stats(w).FunctionTable)
            worker_time = sum([profile_stats(w).FunctionTable.TotalTime]);
            fprintf('Worker %d: Total time = %.2f s\n', w, worker_time);
        end
    end
    fprintf('===========================\n\n');
end
