function varargout = memory_utils(action, varargin)
% Utilities for memory usage tracking across MPI workers
% Tracks maximum memory (RAM) used by each worker at each time step
%
% Usage:
%   memory_utils('init')                          - Initialize memory tracking
%   memory_utils('record', step_num)              - Record memory at current step
%   memory_utils('report', num_workers)           - Collect and display results
%   max_mem = memory_utils('get_max')             - Get maximum memory used (MB)

persistent max_memory_mb;
persistent memory_history;

switch action
    case 'init'
        max_memory_mb = 0;
        memory_history = struct('step', {}, 'memory_mb', {});
        
    case 'record'
        if nargin < 2
            error('memory_utils(''record'') requires step_num argument');
        end
        step_num = varargin{1};
        
        % Get current memory usage
        if ispc
            % Windows: use memory function
            mem_info = memory;
            current_mb = mem_info.MemUsedMATLAB / 1024^2;
        else
            % Unix/Mac: use system command to get actual process memory
            % Get MATLAB's process ID
            pid = feature('getpid');
            % Use ps command to get RSS (Resident Set Size) in KB
            [status, result] = system(sprintf('ps -o rss= -p %d', pid));
            if status == 0
                % RSS is in KB, convert to MB
                current_mb = str2double(strtrim(result)) / 1024;
            else
                % Fallback: use whos (less accurate but better than nothing)
                vars = whos;
                total_bytes = sum([vars.bytes]);
                current_mb = total_bytes / 1024^2;
            end
        end
        
        % Update maximum
        if isempty(max_memory_mb) || current_mb > max_memory_mb
            max_memory_mb = current_mb;
        end
        
        % Store in history
        memory_history(end+1).step = step_num;
        memory_history(end).memory_mb = current_mb;
        
    case 'report'
        if nargin < 2
            error('memory_utils(''report'') requires num_workers argument');
        end
        num_workers = varargin{1};
        
        % Handle case where memory tracking wasn't used
        if isempty(max_memory_mb)
            max_memory_mb = 0;
        end
        
        % Also report workspace memory for comparison
        vars = whos;
        workspace_mb = sum([vars.bytes]) / 1024^2;
        
        fprintf('\n=== Memory Usage Summary ===\n');
        fprintf('Worker %d: Max memory = %.2f MB (workspace: %.2f MB)\n', ...
                spmdIndex, max_memory_mb, workspace_mb);
        
        % Gather max memory from all workers
        max_mem_all = spmdCat(max_memory_mb, 1);
        
        % Only rank 1 prints the full summary
        if spmdIndex == 1
            fprintf('\n--- All Workers ---\n');
            for w = 1:num_workers
                fprintf('Worker %d: Max memory = %.2f MB\n', w, max_mem_all(w));
            end
            fprintf('Total max memory (sum): %.2f MB\n', sum(max_mem_all));
            fprintf('Peak memory per worker (max): %.2f MB\n', max(max_mem_all));
            fprintf('Average memory per worker: %.2f MB\n', mean(max_mem_all));
            fprintf('============================\n\n');
        end
        
    case 'get_max'
        if nargout > 0
            varargout{1} = max_memory_mb;
        end
        
    otherwise
        error('Unknown action: %s', action);
end

end
