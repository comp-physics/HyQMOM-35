function idx = moment_idx(moment_name, order)
%MOMENT_IDX Map moment names to array indices
%   idx = moment_idx('M200')        % Returns index in 35-moment vector
%   idx = moment_idx('M200', 4)     % Explicitly specify order 4
%   idx = moment_idx('S110', 4)     % Works for any moment type (M/C/S)
%   Eliminates magic numbers in code by providing symbolic access:
%     M(moment_idx('M200'))  instead of  M(3)
%   The moment name prefix (M/C/S) is ignored - only indices matter.
%   See also: moment_names, moment_struct
    if length(moment_name) == 3 && all(isstrprop(moment_name, 'digit'))
        moment_name = ['M', moment_name];
    end
    
    % Default to order 4 (35 moments)
    if nargin < 2
        order = 4;
    end
    
    % Get canonical moment names for this order (includes prefix)
    names = moment_names(order);
    
    % Find the index (now comparing with prefix)
    idx = find(strcmp(names, moment_name), 1);
    
    if isempty(idx)
        error('moment_idx:notFound', ...
              'Moment ''%s'' not found in order %d moment list', ...
              moment_name, order);
    end
end

