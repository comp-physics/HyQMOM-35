function varargout = diagnostics(operation, varargin)
%DIAGNOSTICS Unified diagnostic and checking utilities
%
% Syntax:
%   [S30,S40,S11,S21,S31,S12,S22,S03,S13,S04] = diagnostics('check2D', S30,S40,S11,S21,S31,S12,S22,S03,S13,S04)
%   [S300r1, S400r1, ...] = diagnostics('check2D_all_planes', S300, S400, S110, ...)
%   [Diff, MaxDiff] = diagnostics('test_symmetry', M, Np)
%
% Operations:
%   'check2D'           - Check and correct 2D moments in 3D code
%   'check2D_all_planes' - Apply check2D to all three coordinate planes
%   'test_symmetry'     - Check symmetry of moment matrix

    switch operation
        case 'check2D'
            [varargout{1:10}] = check2D_impl(varargin{:});
        case 'check2D_all_planes'
            [varargout{1:30}] = check2D_all_planes_impl(varargin{:});
        case 'test_symmetry'
            [varargout{1:2}] = test_symmetry_impl(varargin{:});
        otherwise
            error('diagnostics:unknownOp', 'Unknown operation: %s', operation);
    end
end

function [S30,S40,S11,S21,S31,S12,S22,S03,S13,S04] = check2D_impl(S30,S40,S11,S21,S31,S12,S22,S03,S13,S04)
% Check and correct 2D moments in 3D code

    h2min = 1000*eps;
    if abs(S11) >= 1
        % Collapse both S11 >= 1 and S11 <= -1 branches
        S11 = sign(S11);
        s3m = sqrt(abs(S30*S03));
        s4m = sqrt(S40*S04);
        H2m = max(h2min, s4m - s3m^2 - 1);
        s4m = H2m + s3m^2 + 1;
        S30 = sign(S30)*s3m;
        S03 = S11*S30;
        S40 = s4m;
        S04 = s4m;
        S12 = S11*S03;
        S21 = S11*S30;
        S13 = S11*S04;
        S31 = S11*S40;
        S22 = s4m;
    else
        % check and correct realizability of cross moments
        [S21,S12,S31,S22,S13] = realizable_2D(S30,S40,S11,S21,S31,S12,S22,S03,S13,S04);
    end
end

function [S300r1, S400r1, S110r, S210r, S310r, S120r, S220r, S030r1, S130r, S040r1, ...
          S300r2, S400r2, S101r, S201r, S301r, S102r, S202r, S003r2, S103r, S004r2, ...
          S030r3, S040r3, S011r, S021r, S031r, S012r, S022r, S003r3, S013r, S004r3] = ...
    check2D_all_planes_impl(S300, S400, S110, S210, S310, S120, S220, S030, S130, S040, ...
                           S101, S201, S301, S102, S202, S003, S103, S004, ...
                           S011, S021, S031, S012, S022, S013)
% Apply check2D to all three coordinate planes

    % Apply check2D to XY plane
    [S300r1, S400r1, S110r, S210r, S310r, S120r, S220r, S030r1, S130r, S040r1] = ...
        check2D_impl(S300, S400, S110, S210, S310, S120, S220, S030, S130, S040);
    
    % Apply check2D to XZ plane
    [S300r2, S400r2, S101r, S201r, S301r, S102r, S202r, S003r2, S103r, S004r2] = ...
        check2D_impl(S300, S400, S101, S201, S301, S102, S202, S003, S103, S004);
    
    % Apply check2D to YZ plane
    [S030r3, S040r3, S011r, S021r, S031r, S012r, S022r, S003r3, S013r, S004r3] = ...
        check2D_impl(S030, S040, S011, S021, S031, S012, S022, S003, S013, S004);
end

function [Diff, MaxDiff] = test_symmetry_impl(M, Np)
% Check symmetry of moment matrix

    Diff = zeros(Np,5);
    for i = 1:Np
        Diff(i,1) = M(i,i,1) - M(Np+1-i,Np+1-i,1);
        Diff(i,2) = M(i,i,2) + M(Np+1-i,Np+1-i,2);
        Diff(i,3) = M(i,i,3) - M(Np+1-i,Np+1-i,3);
        Diff(i,4) = M(i,i,4) + M(Np+1-i,Np+1-i,4);
        Diff(i,5) = M(i,i,5) - M(Np+1-i,Np+1-i,5);
    end
    MaxDiff = zeros(5,1);
    for k = 1:5
        Normk = norm(Diff(:,k));
        MaxDiff(k) = max(Diff(:,k))/(Normk + 1);
    end
end

