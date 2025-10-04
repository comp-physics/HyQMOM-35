function varargout = moment_struct(operation, varargin)
%MOMENT_STRUCT Create and manipulate moment structures
%
%   Provides a clean struct-based interface for moments, eliminating
%   the need for crazy-long argument lists.
%
%   This utility focuses on DATA STRUCTURE management (vector↔struct).
%   For mathematical transformations (S→C, array extraction), use moment_conversion_utils().
%
%   Operations:
%     'from_vector' - Create struct from 35 or 56-element vector
%     'to_vector'   - Convert struct to vector
%     'pack'        - Pack individual moments into struct
%
%   Examples:
%     % Instead of 35 individual variables:
%     M = moment_struct('from_vector', M_vec);  % One struct!
%     
%     % Access individual moments:
%     density = M.M000;
%     vel_x = M.M100 / M.M000;
%     
%     % Pass to functions (much cleaner):
%     [Fx, Fy, Fz] = compute_fluxes(M);
%
%   Complementary Tools:
%     - moment_conversion_utils: Mathematical transformations (S→C, C→M, etc.)
%     - moment_struct: Data structure management (this file)

    switch lower(operation)
        case 'from_vector'
            varargout{1} = vector_to_struct(varargin{:});
        case 'to_vector'
            varargout{1} = struct_to_vector(varargin{:});
        case 'pack'
            varargout{1} = pack_moments(varargin{:});
        case 'unpack'
            % Returns individual variables - use sparingly!
            error('moment_struct:deprecated', ...
                  'Use struct notation instead: M.M000, M.M100, etc.');
        otherwise
            error('moment_struct:UnknownOperation', 'Unknown operation: %s', operation);
    end
end

%% Create struct from vector
function M = vector_to_struct(vec)
%VECTOR_TO_STRUCT Convert moment vector to struct
%
%   Input: 35-element (M4) or 56-element (M5) vector
%   Output: Struct with fields M000, M100, M010, etc.

    % Get canonical moment names from shared utility
    if length(vec) == 35
        names = moment_names(4);
    elseif length(vec) == 56
        names = moment_names(5);
    else
        error('moment_struct:InvalidLength', ...
              'Vector must be 35 (M4) or 56 (M5) elements, got %d', length(vec));
    end
    
    M = struct();
    for i = 1:length(vec)
        M.(names{i}) = vec(i);
    end
end

%% Convert struct to vector
function vec = struct_to_vector(M, order)
%STRUCT_TO_VECTOR Convert moment struct to vector
%
%   Input: 
%     M     - Moment struct
%     order - 4 (35 moments) or 5 (56 moments), default: auto-detect
%   Output: Vector in standard ordering

    if nargin < 2
        % Auto-detect order
        if isfield(M, 'M500')
            order = 5;
        else
            order = 4;
        end
    end
    
    % Get canonical moment names from shared utility
    names = moment_names(order);
    
    vec = zeros(length(names), 1);
    for i = 1:length(names)
        vec(i) = M.(names{i});
    end
end

%% Pack individual moments into struct (for backward compatibility)
function M = pack_moments(varargin)
%PACK_MOMENTS Pack individual moment variables into struct
%
%   Usage: M = moment_struct('pack', 'M000', M000_val, 'M100', M100_val, ...)
%
%   Better: Use from_vector if you have a vector

    M = struct();
    for i = 1:2:length(varargin)
        name = varargin{i};
        value = varargin{i+1};
        M.(name) = value;
    end
end

