function names = moment_names(order)
%MOMENT_NAMES Return canonical moment names for specified order
%
%   names = moment_names(4)  % Returns 35 names for 4th-order moments
%   names = moment_names(5)  % Returns 56 names for 5th-order moments
%
%   This is a shared utility to avoid duplicating moment name lists
%   across moment_conversion_utils.m and moment_struct.m

    % 4th order: 35 moments
    names_35 = {'M000','M100','M200','M300','M400','M010','M110','M210','M310','M020','M120','M220','M030','M130','M040',...
                'M001','M101','M201','M301','M002','M102','M202','M003','M103','M004','M011','M111','M211','M021','M121',...
                'M031','M012','M112','M013','M022'};
    
    % 5th order: 56 moments (35 + 21 closure moments)
    names_21_closure = {'M500','M410','M320','M230','M140','M401','M302','M203','M104','M311','M221','M131','M212','M113','M122',...
                        'M050','M041','M032','M023','M014','M005'};
    
    if order == 4
        names = names_35;
    elseif order == 5
        names = [names_35, names_21_closure];
    else
        error('moment_names:InvalidOrder', 'Order must be 4 or 5, got %d', order);
    end
end

