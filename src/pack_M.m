function M = pack_M(m)
% pack_M Convert struct with named moment fields to moment vector
%
% Input:
%   m - struct with named moment fields (m000, m100, etc.)
%
% Output:
%   M - 35-element moment vector [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,...
%                                 M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,...
%                                 M031,M012,M112,M013,M022]

fields = moment_field_order();
M = zeros(35, 1);

for k = 1:35
    M(k) = m.(fields{k});
end

end
