function [S_out] = check2D_all_planes(S)
% check2D_all_planes Apply check2D to all three coordinate planes
%
% Input:
%   S - standardized moments array (35 elements)
%
% Output:
%   S_out - corrected standardized moments array

S_out = S;

% UV plane (x-y)
[S_out(4),S_out(5),S_out(7),S_out(8),S_out(9),S_out(11),S_out(12),S_out(13),S_out(14),S_out(15)] = ...
    check2D(S(4),S(5),S(7),S(8),S(9),S(11),S(12),S(13),S(14),S(15));

% UW plane (x-z)  
[S_out(4),S_out(5),S_out(17),S_out(18),S_out(19),S_out(21),S_out(22),S_out(23),S_out(24),S_out(25)] = ...
    check2D(S(4),S(5),S(17),S(18),S(19),S(21),S(22),S(23),S(24),S(25));

% VW plane (y-z)
[S_out(13),S_out(15),S_out(26),S_out(29),S_out(31),S_out(32),S_out(35),S_out(23),S_out(34),S_out(25)] = ...
    check2D(S(13),S(15),S(26),S(29),S(31),S(32),S(35),S(23),S(34),S(25));

end
