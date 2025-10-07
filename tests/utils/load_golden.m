function golden_data = load_golden(filename, tests_dir)
%LOAD_GOLDEN Load and validate golden file for testing
%   golden_data = load_golden(filename, tests_dir)
%   
%   Inputs:
%       filename  - Name of golden file (e.g., 'test_closures_golden.mat')
%       tests_dir - Optional path to tests directory (default: auto-detect)
%   
%   Outputs:
%       golden_data - Struct containing golden test data

if nargin < 2
    % Auto-detect tests directory
    current_file = mfilename('fullpath');
    tests_dir = fileparts(fileparts(current_file));
end

% Construct full path to golden file
golden_path = fullfile(tests_dir, 'goldenfiles', filename);

% Check if file exists
if ~exist(golden_path, 'file')
    error('Golden file not found: %s\nRun create_test_goldenfiles.m first.', golden_path);
end

% Load golden data
loaded = load(golden_path);

% Validate structure
if ~isfield(loaded, 'golden_data')
    error('Invalid golden file format: missing golden_data field');
end

golden_data = loaded.golden_data;

% Validate metadata
if ~isfield(golden_data, 'metadata')
    warning('Golden file missing metadata');
else
    fprintf('Loaded golden file: %s\n', filename);
    fprintf('  Version: %s\n', golden_data.metadata.version);
    fprintf('  Created: %s\n', golden_data.metadata.date);
end

% Validate tests field
if ~isfield(golden_data, 'tests')
    error('Invalid golden file format: missing tests field');
end

end
