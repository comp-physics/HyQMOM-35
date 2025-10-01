function [Np, tmax, enable_plots, save_output] = parse_input_args(nargin_val, varargin_cell, defaults)
% PARSE_INPUT_ARGS Parses command-line arguments for main simulation
%
%   [Np, tmax, enable_plots, save_output] = parse_input_args(nargin, varargin, defaults)
%
%   Inputs:
%       nargin_val - Number of input arguments
%       varargin_cell - Cell array of input arguments
%       defaults - Structure with default values
%
%   Outputs:
%       Np - Grid points per dimension
%       tmax - Maximum simulation time
%       enable_plots - Boolean for plotting
%       save_output - Boolean for saving results

switch nargin_val
    case 0
        % Use all defaults
        Np = defaults.Np;
        tmax = defaults.tmax;
        enable_plots = defaults.enable_plots;
        save_output = defaults.save_output;
    case 2
        % Override Np and tmax only
        Np = varargin_cell{1};
        tmax = varargin_cell{2};
        enable_plots = defaults.enable_plots;
        save_output = defaults.save_output;
    case 3
        % Override Np, tmax, and enable_plots
        Np = varargin_cell{1};
        tmax = varargin_cell{2};
        enable_plots = varargin_cell{3};
        save_output = defaults.save_output;
    case 4
        % Override all parameters
        Np = varargin_cell{1};
        tmax = varargin_cell{2};
        enable_plots = varargin_cell{3};
        save_output = varargin_cell{4};
    otherwise
        error(['Invalid number of arguments. Usage:\n' ...
               '  main()\n' ...
               '  main(Np, tmax)\n' ...
               '  main(Np, tmax, enable_plots)\n' ...
               '  main(Np, tmax, enable_plots, save_output)']);
end

end

