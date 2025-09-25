function templates = jacobian_templates(direction)
% jacobian_templates Return index templates for jacobian6 computation
%
% Input:
%   direction - 'x' or 'y' for flux direction
%
% Output:
%   templates - cell array of index templates for moment extraction

switch lower(direction)
    case 'x'
        % X-direction templates
        templates = {
            % UV moments: (m000,m010,m020,m030,m040,m100,m110,m120,m130,m200,m210,m220,m300,m310,m400)
            [1, 6, 10, 13, 15, 2, 7, 11, 14, 3, 8, 12, 4, 9, 5];
            % UW moments: (m000,m001,m002,m003,m004,m100,m101,m102,m103,m200,m201,m202,m300,m301,m400)
            [1, 16, 20, 23, 25, 2, 17, 21, 24, 3, 18, 22, 4, 19, 5]
        };
        
    case 'y'
        % Y-direction templates
        templates = {
            % VU moments: (m000,m100,m200,m300,m400,m010,m110,m210,m310,m020,m120,m220,m030,m130,m040)
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
            % VW moments: (m000,m001,m002,m003,m004,m010,m011,m012,m013,m020,m021,m022,m030,m031,m040)
            [1, 16, 20, 23, 25, 6, 26, 32, 34, 10, 29, 35, 13, 31, 15]
        };
        
    otherwise
        error('jacobian_templates: direction must be ''x'' or ''y''');
end

end
