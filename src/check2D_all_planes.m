function [S300r1, S400r1, S110r, S210r, S310r, S120r, S220r, S030r1, S130r, S040r1, ...
          S300r2, S400r2, S101r, S201r, S301r, S102r, S202r, S003r2, S103r, S004r2, ...
          S030r3, S040r3, S011r, S021r, S031r, S012r, S022r, S003r3, S013r, S004r3] = ...
    check2D_all_planes(S300, S400, S110, S210, S310, S120, S220, S030, S130, S040, ...
                       S101, S201, S301, S102, S202, S003, S103, S004, ...
                       S011, S021, S031, S012, S022, S013)
% check2D_all_planes applies check2D to all three coordinate planes
% This eliminates code duplication from calling check2D three times
%
% Plane 1 (XY): S300, S400, S110, S210, S310, S120, S220, S030, S130, S040
% Plane 2 (XZ): S300, S400, S101, S201, S301, S102, S202, S003, S103, S004
% Plane 3 (YZ): S030, S040, S011, S021, S031, S012, S022, S003, S013, S004

% Apply check2D to XY plane
[S300r1, S400r1, S110r, S210r, S310r, S120r, S220r, S030r1, S130r, S040r1] = ...
    check2D(S300, S400, S110, S210, S310, S120, S220, S030, S130, S040);

% Apply check2D to XZ plane
[S300r2, S400r2, S101r, S201r, S301r, S102r, S202r, S003r2, S103r, S004r2] = ...
    check2D(S300, S400, S101, S201, S301, S102, S202, S003, S103, S004);

% Apply check2D to YZ plane
[S030r3, S040r3, S011r, S021r, S031r, S012r, S022r, S003r3, S013r, S004r3] = ...
    check2D(S030, S040, S011, S021, S031, S012, S022, S003, S013, S004);

end

