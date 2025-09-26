function C5 = scale_S5_to_C5(S5_cell, sC200, sC020, sC002)
% scale_S5_to_C5 Scale standardized 5th-order moments to central moments
%
% Inputs:
%   S5_cell - cell array of 21 standardized 5th-order moments
%   sC200, sC020, sC002 - square roots of diagonal second-order central moments
%
% Output:
%   C5 - cell array of 21 central 5th-order moments

% Exponent table for each moment [x_exp, y_exp, z_exp]
% Order: 500,410,320,230,140,401,302,203,104,311,221,131,212,113,122,050,041,032,023,014,005
exponents = [
    5 0 0; 4 1 0; 3 2 0; 2 3 0; 1 4 0;  % 500-140
    4 0 1; 3 0 2; 2 0 3; 1 0 4;          % 401-104
    3 1 1; 2 2 1; 1 3 1;                 % 311-131
    2 1 2; 1 1 3; 1 2 2;                 % 212-122
    0 5 0; 0 4 1; 0 3 2; 0 2 3; 0 1 4;  % 050-014
    0 0 5                                % 005
];

C5 = cell(size(S5_cell));
for k = 1:numel(S5_cell)
    C5{k} = S5_cell{k} * (sC200^exponents(k,1)) * (sC020^exponents(k,2)) * (sC002^exponents(k,3));
end

end
