function constants = get_constants(Ma)
% get_constants Get simulation constants based on Mach number
%
% Input:
%   Ma - Mach number
%
% Output:
%   constants - struct with simulation constants

constants = struct( ...
    's3max', 4.0 + abs(Ma)/2.0, ...
    'h2min', 1e-8, ...
    'itrealmax', 6 ...
);

end
