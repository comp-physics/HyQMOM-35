function setup_paths(base_dir)
%SETUP_PATHS Add all necessary paths for the HyQMOM simulation
%   setup_paths()            - Auto-detect from current file location
%   setup_paths(base_dir)    - Use specified base directory
%   Adds:
%     - src/
%     - src/autogen/
%   This replaces scattered addpath blocks throughout the codebase.
%   Call once on the client, and once within spmd blocks.
%   Example (in main.m):
%       setup_paths();  % On client
%       spmd
%           setup_paths(base_dir_from_client);  % On workers
%       end
%   See also: addpath
    if nargin < 1
        % Auto-detect base directory from this file's location
        base_dir = fileparts(mfilename('fullpath'));
    end
    
    % Add src directory
    src_dir = fullfile(base_dir, 'src');
    if exist(src_dir, 'dir')
        addpath(src_dir);
        
        % Add autogen subdirectory
        autogen_dir = fullfile(src_dir, 'autogen');
        if exist(autogen_dir, 'dir')
            addpath(autogen_dir);
        end
    else
        warning('setup_paths:srcNotFound', ...
                'Source directory not found: %s', src_dir);
    end
end

