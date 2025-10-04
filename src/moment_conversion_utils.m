function varargout = moment_conversion_utils(operation, varargin)
%MOMENT_CONVERSION_UTILS Utilities for moment array conversions
%
%   Eliminates repetitive moment packing/unpacking code
%
%   Operations:
%     'S_to_C'      : Convert standardized to central moments (vectorized)
%     'M5_to_vars'  : Extract M5 3D array to individual variables
%     'M4_to_vars'  : Extract M4 3D array to individual variables
%
%   Examples:
%     C = moment_conversion_utils('S_to_C', S, sC200, sC020, sC002, indices)
%     [M000, M100, ...] = moment_conversion_utils('M5_to_vars', M5)

    switch lower(operation)
        case 's_to_c'
            varargout{1} = standardized_to_central(varargin{:});
        case 'm5_to_vars'
            [varargout{1:56}] = extract_M5_array(varargin{:});
        case 'm4_to_vars'
            [varargout{1:35}] = extract_M4_array(varargin{:});
        case 'c_to_s'
            varargout{1} = central_to_standardized(varargin{:});
        otherwise
            error('moment_conversion_utils:UnknownOperation', 'Unknown operation: %s', operation);
    end
end

%% Convert standardized to central moments (vectorized)
function C = standardized_to_central(S, sC200, sC020, sC002)
%STANDARDIZED_TO_CENTRAL Vectorized S→C conversion
%
%   C### = S### * sC200^i * sC020^j * sC002^k  where i+j+k = moment order
%
%   Input:
%     S      - Structure or vector of standardized moments
%     sC200  - sqrt(C200)
%     sC020  - sqrt(C020)
%     sC002  - sqrt(C002)
%
%   Output:
%     C - Structure with central moments

% If S is a vector, assume it's ordered like the standard 35-moment vector
if isnumeric(S)
    S_struct = unpack_moments_to_struct(S);
else
    S_struct = S;
end

% Define moment indices: [i, j, k] for C_ijk
% This is the canonical ordering of moments up to 5th order
moment_list = {
    % Order 0-1
    'M000', [0,0,0]; 'M100', [1,0,0]; 'M010', [0,1,0]; 'M001', [0,0,1];
    % Order 2
    'M200', [2,0,0]; 'M110', [1,1,0]; 'M101', [1,0,1];
    'M020', [0,2,0]; 'M011', [0,1,1]; 'M002', [0,0,2];
    % Order 3
    'M300', [3,0,0]; 'M210', [2,1,0]; 'M201', [2,0,1];
    'M120', [1,2,0]; 'M111', [1,1,1]; 'M102', [1,0,2];
    'M030', [0,3,0]; 'M021', [0,2,1]; 'M012', [0,1,2]; 'M003', [0,0,3];
    % Order 4
    'M400', [4,0,0]; 'M310', [3,1,0]; 'M301', [3,0,1];
    'M220', [2,2,0]; 'M211', [2,1,1]; 'M202', [2,0,2];
    'M130', [1,3,0]; 'M121', [1,2,1]; 'M112', [1,1,2]; 'M103', [1,0,3];
    'M040', [0,4,0]; 'M031', [0,3,1]; 'M022', [0,2,2]; 'M013', [0,1,3]; 'M004', [0,0,4];
    % Order 5 (closure moments)
    'M500', [5,0,0]; 'M410', [4,1,0]; 'M320', [3,2,0]; 'M230', [2,3,0]; 'M140', [1,4,0];
    'M401', [4,0,1]; 'M302', [3,0,2]; 'M203', [2,0,3]; 'M104', [1,0,4];
    'M311', [3,1,1]; 'M221', [2,2,1]; 'M131', [1,3,1];
    'M212', [2,1,2]; 'M113', [1,1,3]; 'M122', [1,2,2];
    'M050', [0,5,0]; 'M041', [0,4,1]; 'M032', [0,3,2]; 'M023', [0,2,3]; 'M014', [0,1,4]; 'M005', [0,0,5];
};

C = struct();

% Vectorized conversion
for i = 1:size(moment_list, 1)
    name = moment_list{i,1};
    idx = moment_list{i,2};
    
    % For S→C: C_ijk = S_ijk * sC200^i * sC020^j * sC002^k
    s_name = strrep(name, 'M', 'S');
    if isfield(S_struct, s_name)
        power_factor = (sC200^idx(1)) * (sC020^idx(2)) * (sC002^idx(3));
        C.(strrep(name, 'M', 'C')) = S_struct.(s_name) * power_factor;
    end
end

end

%% Extract M5 array to variables (5th order = 56 moments)
function [M000, M100, M010, M001, M200, M110, M101, M020, M011, M002, ...
          M300, M210, M201, M120, M111, M102, M030, M021, M012, M003, ...
          M400, M310, M301, M220, M211, M202, M130, M121, M112, M103, M040, M031, M022, M013, M004, ...
          M500, M410, M320, M230, M140, M401, M302, M203, M104, M311, M221, M131, M212, M113, M122, ...
          M050, M041, M032, M023, M014, M005] = extract_M5_array(M5)
%EXTRACT_M5_ARRAY Extract 5th-order moment array to individual variables
%
%   Replaces ~60 lines of M### = M5(i,j,k) assignments

% Order 0-1
M000 = M5(1,1,1);
M100 = M5(2,1,1); M010 = M5(1,2,1); M001 = M5(1,1,2);

% Order 2
M200 = M5(3,1,1); M110 = M5(2,2,1); M101 = M5(2,1,2);
M020 = M5(1,3,1); M011 = M5(1,2,2); M002 = M5(1,1,3);

% Order 3
M300 = M5(4,1,1); M210 = M5(3,2,1); M201 = M5(3,1,2);
M120 = M5(2,3,1); M111 = M5(2,2,2); M102 = M5(2,1,3);
M030 = M5(1,4,1); M021 = M5(1,3,2); M012 = M5(1,2,3); M003 = M5(1,1,4);

% Order 4
M400 = M5(5,1,1); M310 = M5(4,2,1); M301 = M5(4,1,2);
M220 = M5(3,3,1); M211 = M5(3,2,2); M202 = M5(3,1,3);
M130 = M5(2,4,1); M121 = M5(2,3,2); M112 = M5(2,2,3); M103 = M5(2,1,4);
M040 = M5(1,5,1); M031 = M5(1,4,2); M022 = M5(1,3,3); M013 = M5(1,2,4); M004 = M5(1,1,5);

% Order 5 (closure)
M500 = M5(6,1,1);
M410 = M5(5,2,1); M320 = M5(4,3,1); M230 = M5(3,4,1); M140 = M5(2,5,1);
M401 = M5(5,1,2); M302 = M5(4,1,3); M203 = M5(3,1,4); M104 = M5(2,1,5);
M311 = M5(4,2,2); M221 = M5(3,3,2); M131 = M5(2,4,2);
M212 = M5(3,2,3); M113 = M5(2,2,4); M122 = M5(2,3,3);
M050 = M5(1,6,1); M041 = M5(1,5,2); M032 = M5(1,4,3); M023 = M5(1,3,4); M014 = M5(1,2,5); M005 = M5(1,1,6);

end

%% Extract M4 array to variables (4th order = 35 moments)
function [M000, M100, M010, M001, M200, M110, M101, M020, M011, M002, ...
          M300, M210, M201, M120, M111, M102, M030, M021, M012, M003, ...
          M400, M310, M301, M220, M211, M202, M130, M121, M112, M103, M040, M031, M022, M013, M004] = extract_M4_array(M4)
%EXTRACT_M4_ARRAY Extract 4th-order moment array to individual variables

% Order 0-1
M000 = M4(1,1,1);
M100 = M4(2,1,1); M010 = M4(1,2,1); M001 = M4(1,1,2);

% Order 2
M200 = M4(3,1,1); M110 = M4(2,2,1); M101 = M4(2,1,2);
M020 = M4(1,3,1); M011 = M4(1,2,2); M002 = M4(1,1,3);

% Order 3
M300 = M4(4,1,1); M210 = M4(3,2,1); M201 = M4(3,1,2);
M120 = M4(2,3,1); M111 = M4(2,2,2); M102 = M4(2,1,3);
M030 = M4(1,4,1); M021 = M4(1,3,2); M012 = M4(1,2,3); M003 = M4(1,1,4);

% Order 4
M400 = M4(5,1,1); M310 = M4(4,2,1); M301 = M4(4,1,2);
M220 = M4(3,3,1); M211 = M4(3,2,2); M202 = M4(3,1,3);
M130 = M4(2,4,1); M121 = M4(2,3,2); M112 = M4(2,2,3); M103 = M4(2,1,4);
M040 = M4(1,5,1); M031 = M4(1,4,2); M022 = M4(1,3,3); M013 = M4(1,2,4); M004 = M4(1,1,5);

end

%% Helper: unpack moments vector to struct (for compatibility)
function S = unpack_moments_to_struct(M_vec)
% Quick unpacker for common 35-moment vector
moment_names = {'S000','S100','S200','S300','S400','S010','S110','S210','S310','S020','S120','S220','S030','S130','S040',...
                'S001','S101','S201','S301','S002','S102','S202','S003','S103','S004','S011','S111','S211','S021','S121',...
                'S031','S012','S112','S013','S022'};
S = struct();
for i = 1:min(length(M_vec), length(moment_names))
    S.(moment_names{i}) = M_vec(i);
end
end

%% Central to standardized (inverse operation)
function S = central_to_standardized(C, sC200, sC020, sC002)
%CENTRAL_TO_STANDARDIZED Vectorized C→S conversion (inverse of S→C)
% For completeness - same logic but inverted

if isnumeric(C)
    C_struct = unpack_moments_to_struct(C);
else
    C_struct = C;
end

moment_list = {
    'M000', [0,0,0]; 'M100', [1,0,0]; 'M010', [0,1,0]; 'M001', [0,0,1];
    'M200', [2,0,0]; 'M110', [1,1,0]; 'M101', [1,0,1]; 'M020', [0,2,0]; 'M011', [0,1,1]; 'M002', [0,0,2];
    'M300', [3,0,0]; 'M210', [2,1,0]; 'M201', [2,0,1]; 'M120', [1,2,0]; 'M111', [1,1,1]; 'M102', [1,0,2];
    'M030', [0,3,0]; 'M021', [0,2,1]; 'M012', [0,1,2]; 'M003', [0,0,3];
    'M400', [4,0,0]; 'M310', [3,1,0]; 'M301', [3,0,1]; 'M220', [2,2,0]; 'M211', [2,1,1]; 'M202', [2,0,2];
    'M130', [1,3,0]; 'M121', [1,2,1]; 'M112', [1,1,2]; 'M103', [1,0,3];
    'M040', [0,4,0]; 'M031', [0,3,1]; 'M022', [0,2,2]; 'M013', [0,1,3]; 'M004', [0,0,4];
};

S = struct();
for i = 1:size(moment_list, 1)
    name = moment_list{i,1};
    idx = moment_list{i,2};
    c_name = strrep(name, 'M', 'C');
    if isfield(C_struct, c_name)
        power_factor = (sC200^idx(1)) * (sC020^idx(2)) * (sC002^idx(3));
        S.(strrep(name, 'M', 'S')) = C_struct.(c_name) / power_factor;
    end
end

end

