function varargout = realizability(operation, varargin)
%REALIZABILITY Consolidated utility functions for moment realizability checks.
%   This function acts as a dispatcher for various realizability helpers.
%
%   Usage:
%     [S21, S12, S31, S22, S13] = realizability('2D', S30, S40, S11, S21, S31, S12, S22, S03, S13, S04)
%     [S300,...,flag220] = realizability('3D', S300, ..., S022)
%     [S111r] = realizability('S111', S110, S101, S011, S210, S201, S120, S021, S102, S012, S111)
%     [S110r, S101r, S011r, S2r] = realizability('S2', S110, S101, S011)
%     [S210r, S201r] = realizability('S210', S110, S101, S011, S300, S210, S201, H200, beta)
%     S211r = realizability('S211', e11, e22, e33, e12, e13, d23, S211, beta)
%     [S310r, S220r] = realizability('S310', S110, S101, S011, S300, S210, S201, S120, S111, S310, S220, H200, beta)
%     [S220r] = realizability('S310_220', S110, S101, S011, S210, S120, S111, S220)
%     S220r = realizability('S220', S110, S220, A220)
%
%   Operations:
%     '2D'       : Find realizable moment set in 2D
%     '3D'       : Check and correct realizability of cross moments in 3D
%     'S111'     : Check and correct realizability of S111
%     'S2'       : Check and correct realizability of 2nd-order moments
%     'S210'     : Check and correct realizability of S210, S201
%     'S211'     : Check and correct realizability of S211
%     'S310'     : Check and correct realizability of S310 and S220
%     'S310_220' : Check and correct realizability of S220 for S310
%     'S220'     : Check maximum bounds and correct S220
%
%   See also: Flux_closure35_and_realizable_3D

    switch lower(operation)
        case '2d'
            [varargout{1:5}] = realizable_2D(varargin{:});
        case '3d'
            [varargout{1:29}] = realizable_3D(varargin{:});
        case 's111'
            varargout{1} = realizability_S111(varargin{:});
        case 's2'
            [varargout{1:4}] = realizability_S2(varargin{:});
        case 's210'
            [varargout{1:2}] = realizability_S210(varargin{:});
        case 's211'
            varargout{1} = realizability_S211(varargin{:});
        case 's310'
            [varargout{1:2}] = realizability_S310(varargin{:});
        case 's310_220'
            varargout{1} = realizability_S310_220(varargin{:});
        case 's220'
            varargout{1} = realizablity_S220(varargin{:});
        otherwise
            error('realizability:UnknownOperation', 'Unknown operation: %s', operation);
    end
end

%% realizable_2D
function [S21,S12,S31,S22,S13] = realizable_2D(S30,S40,S11,S21,S31,S12,S22,S03,S13,S04)
%realizable_2D Find realizable moment set in 2D
%   
% check realizability of S11
Del1= max(0,1 - S11^2);
%
H20 = max(eps,S40 - S30^2 - 1);
H02 = max(eps,S04 - S03^2 - 1);
% check realizability of S12 and S21
G1 = sqrt(Del1*H02);
s12min = S11*S03 - G1;
s12max = S11*S03 + G1;
if S12 <= s12min 
    S12 = s12min ;
elseif S12 >= s12max
    S12 = s12max;
end
%
G1 = sqrt(Del1*H20);
s21min = S11*S30 - G1;
s21max = S11*S30 + G1;
if S21 <= s21min
    S21 = s21min;
elseif S21 >= s21max
    S21 = s21max;
end
% at this point first minor is nonnegative

% check realizability of S22
G22 = sqrt((H20+S30^2)*(H02+S03^2));
s22min = max(S11^2,1-G22);
% s22min = max(s22min,S11^2 + (S21^2 - 2*S11*S12*S21 + S12^2)/(Del1+eps));
s22max = 1+G22;
if S22 < s22min
    S22 = s22min;
elseif S22 > s22max
    S22 = s22max;
end

% given s22, check realizability of S13 and S31
G31 = (Del1*S22 - Del1*S11^2 - S12^2 + 2*S11*S12*S21 - S21^2)*(Del1*H20 - (S21 - S11*S30)^2);
G13 = (Del1*S22 - Del1*S11^2 - S12^2 + 2*S11*S12*S21 - S21^2)*(Del1*H02 - (S12 - S11*S03)^2);
if G31 < 0 || G13 < 0
    G31 = 0;
    G13 = 0;
else
    G13 = sqrt(G13);
    G31 = sqrt(G31);
end

s31min = S11 + (S12*S21 + S21*S30 - S11*S21^2 - S11*S12*S30 - G31)/(Del1+eps);
s31max = S11 + (S12*S21 + S21*S30 - S11*S21^2 - S11*S12*S30 + G31)/(Del1+eps);
if S31 <= s31min
    S31 = s31min;
elseif S31 >= s31max
    S31 = s31max;
end
 
s13min = S11 + (S12*S21 + S12*S03 - S11*S12^2 - S11*S21*S03 - G13)/(Del1+eps);
s13max = S11 + (S12*S21 + S12*S03 - S11*S12^2 - S11*S21*S03 + G13)/(Del1+eps);
if S13 <= s13min
    S13 = s13min;
elseif S13 >= s13max
    S13 = s13max;
end

% at this point second minor is nonnegative, check third for S22
L3 = delta2starchol_L3(S03,S04,S11,S12,S13,S21,S22,S30,S31,S40);
if L3 < 0
    % check realizability of S22 using roots of degree 3 polynomial
    R = rootsR(Del1,H02,H20,S03,S11,S12,S13,S21,S30,S31);
    Rr = sort(real(R));
    s22maxr = min(s22max,Rr(3));
    s22minr = max(s22min,Rr(2));
    if s22maxr < s22minr
        s22minr = s22maxr;
    end
    S22a = S22;
    if S22 > s22maxr
        S22a = s22maxr;
    elseif S22 < s22minr
        S22a = s22minr;
    end
    % check for complex roots, in which case L3 remains < 0
    if max(abs(imag(R)))/max(abs(R)) > 1000*eps
        % complex roots : do not change S22
    else
        S22 = S22a;
    end
end
end

%% realizable_3D
function [S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
          S101,S201,S301,S102,S202,S003,S103,S004,S011,S111,...
          S211,S021,S121,S031,S012,S112,S013,S022,flag220] = ...
   realizable_3D(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
                 S101,S201,S301,S102,S202,S003,S103,S004,S011,S111,...
                 S211,S021,S121,S031,S012,S112,S013,S022)
% realizable_3D checks and correct realizability of cross moments
 
flag220 = 0;
diagmin = 1.d-10;
h2min = 1.d-10;
S2min = 1.d-12;
 
H200 = max([eps S400-S300^2-1]);
H020 = max([eps S040-S030^2-1]);
H002 = max([eps S004-S003^2-1]);
 
% 4-order moments: check maximum bounds on S220, S202, S022
A220 = sqrt((H200+S300^2)*(H020+S030^2));
S220max = realizablity_S220(S110,S220,A220);
A202 = sqrt((H200+S300^2)*(H002+S003^2));
S202max = realizablity_S220(S101,S202,A202);
A022 = sqrt((H020+S030^2)*(H002+S003^2));
S022max = realizablity_S220(S011,S022,A022);
S220 = min(S220,S220max);
S202 = min(S202,S202max);
S022 = min(S022,S022max);
 
% check and correct realizability of S110, S101, S011
[S110,S101,S011,S2] = realizability_S2(S110,S101,S011);
% S2 = max(0,S2);
% these R should all be positive!
R110 = max(0,1-S110^2);
R110 = sqrt(R110);
R101 = max(0,1-S101^2);
R101 = sqrt(R101);
R011 = max(0,1-S011^2);
R011 = sqrt(R011);
if R110*R101*R011 == 0
    warning('R110*R101*R011 == 0')
end
% S2 = max(0,1 + 2*S110*S101*S011 - S110^2 - S101^2 - S011^2);
if S2 <= S2min % treat faces
    disp('treat faces')
    Rmax = max([R110 R101 R011]);
    if R110 == Rmax  % Sij0 is known
        gam1 = (S101-S110*S011)/(1-S110^2);
        gam2 = (S011-S110*S101)/(1-S110^2);
        S021= gam2*S030 + gam1*S120;
        S201= gam2*S210 + gam1*S300;
        S031= gam2*S040 + gam1*S130;
        S301= gam2*S310 + gam1*S400;
        S111= gam2*S120 + gam1*S210;
        S211= gam2*S220 + gam1*S310;
        S121= gam2*S130 + gam1*S220;
        %
        S012= gam2^2*S030 + 2*gam1*gam2*S120 + gam1^2*S210;
        S102= gam2^2*S120 + 2*gam1*gam2*S210 + gam1^2*S300;
        S022= gam2^2*S040 + 2*gam1*gam2*S130 + gam1^2*S220;
        S202= gam2^2*S220 + 2*gam1*gam2*S310 + gam1^2*S400;
        S112= gam2^2*S130 + 2*gam1*gam2*S220 + gam1^2*S310;
        %
        S003= gam2^3*S030 + 3*gam1*gam2^2*S120 + 3*gam1^2*gam2*S210 + gam1^3*S300;
        S013= gam2^3*S040 + 3*gam1*gam2^2*S130 + 3*gam1^2*gam2*S220 + gam1^3*S310;
        S103= gam2^3*S130 + 3*gam1*gam2^2*S220 + 3*gam1^2*gam2*S310 + gam1^3*S400;
        %
        S004= gam2^4*S040 + 4*gam1*gam2^3*S130 + 6*gam1^2*gam2^2*S220 + 4*gam1^3*gam2*S310 + gam1^4*S400;
    elseif R101 == Rmax   % Si0k is known
        gam3 = (S110-S101*S011)/(1-S101^2);
        gam4 = (S011-S110*S101)/(1-S101^2);
        S012= gam4*S003 + gam3*S102;
        S210= gam4*S201 + gam3*S300;
        S013= gam4*S004 + gam3*S103;
        S310= gam4*S301 + gam3*S400;
        S111= gam4*S102 + gam3*S201;
        S211= gam4*S202 + gam3*S301;
        S112= gam4*S103 + gam3*S202;
        %
        S021= gam4^2*S003 + 2*gam3*gam4*S102 + gam3^2*S201;
        S120= gam4^2*S102 + 2*gam3*gam4*S201 + gam3^2*S300;
        S022= gam4^2*S004 + 2*gam3*gam4*S103 + gam3^2*S202;
        S220= gam4^2*S202 + 2*gam3*gam4*S301 + gam3^2*S400;
        S121= gam4^2*S103 + 2*gam3*gam4*S202 + gam3^2*S301;
        %
        S030= gam4^3*S003 + 3*gam3*gam4^2*S102 + 3*gam3^2*gam4*S201 + gam3^3*S300;
        S031= gam4^3*S004 + 3*gam3*gam4^2*S103 + 3*gam3^2*gam4*S202 + gam3^3*S301;
        S130= gam4^3*S103 + 3*gam3*gam4^2*S202 + 3*gam3^2*gam4*S301 + gam3^3*S400;
        %
        S040= gam4^4*S004 + 4*gam3*gam4^3*S103 + 6*gam3^2*gam4^2*S202 + 4*gam3^3*gam4*S301 + gam3^4*S400;
    else    % S0jk is known
        gam5 = (S101-S110*S011)/(1-S011^2);
        gam6 = (S110-S101*S011)/(1-S011^2);
        S102= gam6*S003 + gam5*S012;
        S120= gam6*S021 + gam5*S030;
        S103= gam6*S004 + gam5*S013;
        S130= gam6*S031 + gam5*S040;
        S111= gam6*S012 + gam5*S021;
        S121= gam6*S022 + gam5*S031;
        S112= gam6*S013 + gam5*S022;
        %
        S201= gam6^2*S003 + 2*gam5*gam6*S012 + gam5^2*S021;
        S210= gam6^2*S012 + 2*gam5*gam6*S021 + gam5^2*S030;
        S202= gam6^2*S004 + 2*gam5*gam6*S013 + gam5^2*S022;
        S220= gam6^2*S022 + 2*gam5*gam6*S031 + gam5^2*S040;
        S211= gam6^2*S013 + 2*gam5*gam6*S022 + gam5^2*S031;
        %
        S300= gam6^3*S003 + 3*gam5*gam6^2*S012 + 3*gam5^2*gam6*S021 + gam5^3*S030;
        S301= gam6^3*S004 + 3*gam5*gam6^2*S013 + 3*gam5^2*gam6*S022 + gam5^3*S031;
        S310= gam6^3*S013 + 3*gam5*gam6^2*S022 + 3*gam5^2*gam6*S031 + gam5^3*S040;
        %
        S400= gam6^4*S004 + 4*gam5*gam6^3*S013 + 6*gam5^2*gam6^2*S022 + 4*gam5^3*gam6*S031 + gam5^4*S040;
    end
     
else % treat interior of 2nd-order moment space
    beta = 1;
     
    %% check for for too small diagonal elements
    E1 = delta2star3D(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
                      S101,S201,S301,S102,S202,S003,S103,S004,S011,S111,...
                      S211,S021,S121,S031,S012,S112,S013,S022);
     
    % 3-order moments: first check and correct realizability of S210, etc.
    flagE11 = 0;
    if E1(1,1) < 0
        [S210,S201] = realizability_S210(S110,S101,S011,S300,S210,S201,H200,beta);
        flagE11 = 1;
    end
    % check and correct realizability of S120, S021
    flagE44 = 0;
    if E1(4,4) < 0
        [S120,S021] = realizability_S210(S110,S011,S101,S030,S120,S021,H020,beta);
        flagE44 = 1;
    end
    % check and correct realizability of S102, S012
    flagE66 = 0;
    if E1(6,6) < 0
        [S012,S102] = realizability_S210(S011,S101,S110,S003,S012,S102,H002,beta);
        flagE66 = 1;
    end
    % check and correct realizability of S111
    [S111r] = realizability_S111(S110,S101,S011,S210,S201,S120,S021,S102,S012,S111);
    if S111r ~= S111
        S111 = S111r;
    end
    % set minimum value based on S310
    S220_310 = realizability_S310_220(S110,S101,S011,S210,S120,S111,S220);
    S202_310 = realizability_S310_220(S101,S110,S011,S201,S102,S111,S202);
    S022_310 = realizability_S310_220(S011,S110,S101,S021,S012,S111,S022);
    %
    % 4-order moments: lower bound for positive diagonal elements
    S220_diag = S220;
    S202_diag = S202;
    S022_diag = S022;
     
    S22min = lower_bound_S220(S011,S012,S021,S101,S102,S110,S111,S120,S201,S210);
     
    flagE22 = 0;
    if S220 < S22min(1)
        S220_diag = S22min(1);
        flagE22 = 1;
        flag220 = 1;
    end
    flagE33 = 0;
    if S202 < S22min(2)
        S202_diag = S22min(2);
        flagE33 = 1;
        flag220 = 1;
    end
    flagE55 = 0;
    if S022 < S22min(3)
        S022_diag = S22min(3);
        flagE55 = 1;
        flag220 = 1;
    end
     
    % recheck for for negative diagonal elements
    E1 = delta2star3D(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
                      S101,S201,S301,S102,S202,S003,S103,S004,S011,S111,...
                      S211,S021,S121,S031,S012,S112,S013,S022);
     
    % check for too small diagonal elements
    if E1(1,1) < diagmin
        flagE11 = 1;
    end
    if E1(2,2) < diagmin
        flagE22 = 1;
    end
    if E1(3,3) < diagmin
        flagE33 = 1;
    end
    if E1(4,4) < diagmin
        flagE44 = 1;
    end
    if E1(5,5) < diagmin
        flagE55 = 1;
    end
    if E1(6,6) < diagmin
        flagE66 = 1;
    end
    
    % compute bound values for non-diagonal moments
    if flagE11 == 1 || flagE22 == 1 || flagE33 == 1 || flagE44 == 1 || flagE55 == 1 || flagE66 == 1
        Mbound = bound_minor1(S003,S011,S012,S021,S030,S101,S102,S110,S111,S120,S201,S210,S300);

        % 4th-order cross moments to 3rd-order boundary values
        S310b = Mbound(1);
        S301b = Mbound(2);
        S130b = Mbound(3);
        S103b = Mbound(4); 
        S031b = Mbound(5); 
        S013b = Mbound(6);
        S211b = Mbound(10);
        S121b = Mbound(11); 
        S112b = Mbound(12);

        % treat cases where a diagonal element is near of below diagmin
        if flagE11 == 1
            S310 = S310b;
            S301 = S301b;
            S211 = S211b;
        end
        if flagE22 == 1
            S310 = S310b;
            S130 = S130b;
            S211 = S211b;
            S121 = S121b;
            S112 = S112b;
        end
        if flagE33 == 1
            S301 = S301b;
            S103 = S103b;
            S211 = S211b;
            S121 = S121b;
            S112 = S112b;
        end
        if flagE44 == 1
            S130 = S130b;
            S031 = S031b;
            S121 = S121b;
        end
        if flagE55 == 1
            S031 = S031b;
            S013 = S013b;
            S211 = S211b;
            S121 = S121b;
            S112 = S112b;
        end
        if flagE66 == 1
            S103 = S103b;
            S013 = S013b;
            S112 = S112b;
        end
    end

    % return if all off-diagonal terms have been fixed
    if flagE11 == 1 && flagE22 == 1 && flagE33 == 1 && flagE44 == 1 && flagE55 == 1 && flagE66 == 1
        disp('all flags are 1')
        return
    end
    
    % Move to fourth-order moments
    % compute varDelta_2^* matrices
    [E1,E2,E3,E4,E5,E6] = delta2star3D_permutation(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
       S101,S201,S301,S102,S202,S003,S103,S004,S011,S111,S211,S021,S121,S031,S012,S112,S013,S022);
     
    % check and correct realizability of S310, etc. if they are still free
    flag310 = 0;
    S220a = S220;
    S220b = S220;
    S202a = S202;
    S202b = S202;
    S022a = S022;
    S022b = S022;
    %
    beta = 1-1000*eps;
    if flagE11 == 0 && flagE22 == 0
        if det(E1(1:2,1:2)) < h2min
            [S310,S220a]= realizability_S310(S110,S101,S011,S300,S210,S201,S120,S111,S310,S220,H200,beta);
            flag310 = 1;
        end
    end
    if flagE22 == 0 && flagE44 == 0
        if det(E3(1:2,1:2)) < h2min
            [S130,S220b]= realizability_S310(S110,S011,S101,S030,S120,S021,S210,S111,S130,S220,H020,beta);
            flag310 = 1;
        end
    end
    if flagE11 == 0 && flagE33 == 0
        if det(E2(1:2,1:2)) < h2min
            [S301,S202a]= realizability_S310(S101,S110,S011,S300,S201,S210,S102,S111,S301,S202,H200,beta);
            flag310 = 1;
        end
    end
    if flagE33 == 0 && flagE66 == 0
        if det(E4(1:2,1:2)) < h2min
            [S103,S202b]= realizability_S310(S101,S011,S110,S003,S102,S012,S201,S111,S103,S202,H002,beta);
            flag310 = 1;
        end
    end
    if flagE44 == 0 && flagE55 == 0
        if det(E5(1:2,1:2)) < h2min
            [S031,S022a]= realizability_S310(S011,S110,S101,S030,S021,S120,S012,S111,S031,S022,H020,beta);
            flag310 = 1;
        end
    end
    if flagE55 == 0 && flagE66 == 0
        if det(E6(1:2,1:2)) < h2min
            [S013,S022b]= realizability_S310(S011,S101,S110,S003,S012,S102,S021,S111,S013,S022,H002,beta);
            flag310 = 1;
        end
    end
    %
    S220r = max(S220a,S220b);
    S220_310r = S220_310;
    if S220 < S220r
        S220_310r = S220r;
        flag220 = 1;
    end
    S022r = max(S022a,S022b);
    S022_310r = S022_310;
    if S022 < S022r
        S022_310r = S022r;
        flag220 = 1;
    end
    S202r = max(S202a,S202b);
    S202_310r = S202_310;
    if S202 < S202r
        S202_310r = S202r;
        flag220 = 1;
    end
     
    % if moments have been corrected, recompute matrices
    if flag310 == 1
        [E1,E2,E3,E4,E5,E6] = delta2star3D_permutation(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
           S101,S201,S301,S102,S202,S003,S103,S004,S011,S111,S211,S021,S121,S031,S012,S112,S013,S022);
        if min([det(E1(1:2,1:2)) det(E2(1:2,1:2)) det(E3(1:2,1:2)) det(E4(1:2,1:2)) det(E5(1:2,1:2)) det(E6(1:2,1:2))]) < 0
            flag220 = 1;
        end
    end
    
    S220a=S220;
    S220b=S220;
    S202a=S202;
    S202c=S202;
    S022b=S022;
    S022c=S022;

    % check and correct realizability of S211, etc. if they are still free
    if flagE11 == 0 && flagE22 == 0 && flagE33 == 0 && flagE55 == 0  
        if det(E1(1:3,1:3)) <= h2min || det(E2(1:3,1:3)) <= h2min
            e11a = E1(1,1);
            e12a = E1(1,2);
            e13a = E1(1,3);
            e22a = E1(2,2);
            e23a = E1(2,3);
            e33a = E1(3,3);
            d23a = S211-e23a;
            d22a = S220-e22a;
            d33a = S202-e33a;
            if e11a > diagmin
                S220a = d22a+e13a^2/e11a;
                S202a = d33a+e12a^2/e11a;
            end
            S211 = realizability_S211(e11a,e22a,e33a,e12a,e13a,d23a,S211,beta);
            flag220 = 1;
        end
    end
    if flagE22 == 0 && flagE33 == 0 && flagE44 == 0 && flagE55 == 0 
        if det(E3(1:3,1:3)) < h2min || det(E5(1:3,1:3)) < h2min
            e11b = E3(1,1);
            e12b = E3(1,2);
            e13b = E3(1,3);
            e22b = E3(2,2);
            e23b = E3(2,3);
            e33b = E3(3,3);
            d23b = S121-e23b;
            d22b = S220-e22b;
            d33b = S022-e33b;
            if e11b > diagmin
                S220b = d22b+e13b^2/e11b;
                S022b = d33b+e12b^2/e11b;
            end
            S121 = realizability_S211(e11b,e22b,e33b,e12b,e13b,d23b,S121,beta);
            flag220 = 1;
        end
    end
    if flagE22 == 0 && flagE33 == 0 && flagE55 == 0 && flagE66 == 0
        % check and correct realizability of S112 given S202 and S022
        if det(E4(1:3,1:3)) < h2min || det(E6(1:3,1:3)) < h2min
            e11c = E6(1,1);
            e12c = E6(1,2);
            e13c = E6(1,3);
            e22c = E6(2,2);
            e23c = E6(2,3);
            e33c = E6(3,3);
            d23c = S112-e23c;
            d22c = S022-e22c;
            d33c = S202-e33c;
            if e11c > diagmin
                S202c = d22c+e13c^2/e11c;
                S022c = d33c+e12c^2/e11c;
            end
            %
            S112 = realizability_S211(e11c,e22c,e33c,e12c,e13c,d23c,S112,beta);
            %
            flag220 = 1;
        end
    end
     
    S220_211 = max(S220a,S220b);
    S202_211 = max(S202a,S202c);
    S022_211 = max(S022b,S022c);
     
    %% CONTINUE FOR 4TH MINOR to correct S220, etc.
     
    E1 = delta2star3D(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
                      S101,S201,S301,S102,S202,S003,S103,S004,S011,S111,...
                      S211,S021,S121,S031,S012,S112,S013,S022);
     
    E3 = delta2star3D(S030,S040,S110,S120,S130,S210,S220,S300,S310,S400,...
                      S011,S021,S031,S012,S022,S003,S013,S004,S101,S111,...
                      S121,S201,S211,S301,S102,S112,S103,S202);
     
    s220mina = S220;
    s220maxa = S220;
    if det(E1(1:4,1:4)) < 0
        flag220 = 1;
        E = E1;
        %
        e11=E(1,1);
        e12=E(1,2);
        e13=E(1,3);
        e14=E(1,4);
        e15=E(1,5);
        e16=E(1,6);
        e22=E(2,2);
        e23=E(2,3);
        e24=E(2,4);
        e25=E(2,5);
        e26=E(2,6);
        e33=E(3,3);
        e34=E(3,4);
        e35=E(3,5);
        e36=E(3,6);
        e44=E(4,4);
        e45=E(4,5);
        e46=E(4,6);
        e55=E(5,5);
        e56=E(5,6);
        e66=E(6,6);
        %
        ex = -e22+e14;
        d22 = S220-e22;
        %
        Y = e33;
        Ra = rootsR_X_Y(Y,e11,e12,e13,e23,e24,e34,e44,ex);
        Rr = sort(real(Ra));
        if max(abs(imag(Ra)))/max(abs(Ra)) < 1000*eps
            s220mina = Rr(2)+d22;
            s220maxa = Rr(3)+d22;
        end
    end
    s220minb = S220;
    s220maxb = S220;
    if det(E3(1:4,1:4)) < 0
        flag220 = 1;
        E = E3;
        %
        e11=E(1,1);
        e12=E(1,2);
        e13=E(1,3);
        e14=E(1,4);
        e15=E(1,5);
        e16=E(1,6);
        e22=E(2,2);
        e23=E(2,3);
        e24=E(2,4);
        e25=E(2,5);
        e26=E(2,6);
        e33=E(3,3);
        e34=E(3,4);
        e35=E(3,5);
        e36=E(3,6);
        e44=E(4,4);
        e45=E(4,5);
        e46=E(4,6);
        e55=E(5,5);
        e56=E(5,6);
        e66=E(6,6);
         
        ex = -e22+e14;
        d22 = S220-e22;
         
        Y = e33;
        Ra = rootsR_X_Y(Y,e11,e12,e13,e23,e24,e34,e44,ex);
        Rr = sort(real(Ra));
        if max(abs(imag(Ra)))/max(abs(Ra)) < 1000*eps
            s220minb = Rr(2)+d22;
            s220maxb = Rr(3)+d22;
        end
    end
     
    E2 = delta2star3D(S300,S400,S101,S201,S301,S102,S202,S003,S103,S004,...
                  S110,S210,S310,S120,S220,S030,S130,S040,S011,S111,...
                  S211,S012,S112,S013,S021,S121,S031,S022);
     
    E4 = delta2star3D(S003,S004,S101,S102,S103,S201,S202,S300,S301,S400,...
                  S011,S012,S013,S021,S022,S030,S031,S040,S110,S111,...
                  S112,S210,S211,S310,S120,S121,S130,S220);
     
    s202mina = S202;
    s202maxa = S202;
    if det(E2(1:4,1:4)) < 0
        flag220 = 1;
        E = E2;
         
        e11=E(1,1);
        e12=E(1,2);
        e13=E(1,3);
        e14=E(1,4);
        e15=E(1,5);
        e16=E(1,6);
        e22=E(2,2);
        e23=E(2,3);
        e24=E(2,4);
        e25=E(2,5);
        e26=E(2,6);
        e33=E(3,3);
        e34=E(3,4);
        e35=E(3,5);
        e36=E(3,6);
        e44=E(4,4);
        e45=E(4,5);
        e46=E(4,6);
        e55=E(5,5);
        e56=E(5,6);
        e66=E(6,6);
         
        ex = -e22+e14;
        d22 = S202-e22;
         
        Y = e33;
        Ra = rootsR_X_Y(Y,e11,e12,e13,e23,e24,e34,e44,ex);
        Rr = sort(real(Ra));
        if max(abs(imag(Ra)))/max(abs(Ra)) < 1000*eps
            s202mina = Rr(2)+d22;
            s202maxa = Rr(3)+d22;
        end
    end
    s202minb = S202;
    s202maxb = S202;
    if det(E4(1:4,1:4)) < 0
        flag220 = 1;
        E = E4;
        %
        e11=E(1,1);
        e12=E(1,2);
        e13=E(1,3);
        e14=E(1,4);
        e15=E(1,5);
        e16=E(1,6);
        e22=E(2,2);
        e23=E(2,3);
        e24=E(2,4);
        e25=E(2,5);
        e26=E(2,6);
        e33=E(3,3);
        e34=E(3,4);
        e35=E(3,5);
        e36=E(3,6);
        e44=E(4,4);
        e45=E(4,5);
        e46=E(4,6);
        e55=E(5,5);
        e56=E(5,6);
        e66=E(6,66);
         
        ex = -e22+e14;
        d22 = S202-e22;
         
        Y = e33;
        Ra = rootsR_X_Y(Y,e11,e12,e13,e23,e24,e34,e44,ex);
        Rr = sort(real(Ra));
        if max(abs(imag(Ra)))/max(abs(Ra)) < 1000*eps
            s202minb = Rr(2)+d22;
            s202maxb = Rr(3)+d22;
        end   
    end
    %
    E5 = delta2star3D(S030,S040,S011,S021,S031,S012,S022,S003,S013,S004,...
                  S110,S120,S130,S210,S220,S300,S310,S400,S101,S111,...
                  S121,S102,S112,S103,S201,S211,S301,S202);
    %
    E6 = delta2star3D(S003,S004,S011,S012,S013,S021,S022,S030,S031,S040,...
                  S101,S102,S103,S201,S202,S300,S301,S400,S110,S111,...
                  S112,S120,S121,S130,S210,S211,S310,S220);
    %
    s022mina = S022;
    s022maxa = S022;
    if det(E5(1:4,1:4)) < 0
        flag220 = 1;
        E = E5;
        %
        e11=E(1,1);
        e12=E(1,2);
        e13=E(1,3);
        e14=E(1,4);
        e15=E(1,5);
        e16=E(1,6);
        e22=E(2,2);
        e23=E(2,3);
        e24=E(2,4);
        e25=E(2,5);
        e26=E(2,6);
        e33=E(3,3);
        e34=E(3,4);
        e35=E(3,5);
        e36=E(3,6);
        e44=E(4,4);
        e45=E(4,5);
        e46=E(4,6);
        e55=E(5,5);
        e56=E(5,6);
        e66=E(6,6);
        %
        ex = -e22+e14;
        d22 = S022-e22;
        %
        Y = e33;
        Ra = rootsR_X_Y(Y,e11,e12,e13,e23,e24,e34,e44,ex);
        Rr = sort(real(Ra));
        if max(abs(imag(Ra)))/max(abs(Ra)) < 1000*eps
            s022mina = Rr(2)+d22;
            s022maxa = Rr(3)+d22;
        end
    end
    s022minb = S022;
    s022maxb = S022;
    if det(E6(1:4,1:4)) < 0
        flag220 = 1;
        E = E6;
        %
        e11=E(1,1);
        e12=E(1,2);
        e13=E(1,3);
        e14=E(1,4);
        e15=E(1,5);
        e16=E(1,6);
        e22=E(2,2);
        e23=E(2,3);
        e24=E(2,4);
        e25=E(2,5);
        e26=E(2,6);
        e33=E(3,3);
        e34=E(3,4);
        e35=E(3,5);
        e36=E(3,6);
        e44=E(4,4);
        e45=E(4,5);
        e46=E(4,6);
        e55=E(5,5);
        e56=E(5,6);
        e66=E(6,6);
        %
        ex = -e22+e14;
        d22 = S022-e22;
        %
        Y = e33;
        Ra = rootsR_X_Y(Y,e11,e12,e13,e23,e24,e34,e44,ex);
        Rr = sort(real(Ra));
        if max(abs(imag(Ra)))/max(abs(Ra)) < 1000*eps
            s022minb = Rr(2)+d22;
            s022maxb = Rr(3)+d22;
        end
    end
    %
    S220 = max([s220mina,s220minb,S220_diag,S220_310,S220_310r,S220_211]);
    S220 = min([S220,S220max,s220maxa,s220maxb]);
    S202 = max([s202mina,s202minb,S202_diag,S202_310,S202_310r,S202_211]);
    S202 = min([S202,S202max,s202maxa,s202maxb]);
    S022 = max([s022mina,s022minb,S022_diag,S022_310,S022_310r,S022_211]);
    S022 = min([S022,S022max,s022maxa,s022maxb]);
end
%
end

%% realizability_S111
function [S111r] = realizability_S111(S110,S101,S011,S210,S201,S120,S021,S102,S012,S111)
% realizability_S111 checks and corrects realizablity of S111
%
S111r = S111;
A110 = ((S101-S011*S110)*S210+(S011-S101*S110)*S120)/(1-S110^2);
A101 = ((S110-S011*S101)*S201+(S011-S110*S101)*S102)/(1-S101^2);
A011 = ((S110-S101*S011)*S021+(S101-S110*S011)*S012)/(1-S011^2);
Rmin = min([A110,A101,A011]);
Rmax = max([A110,A101,A011]);
if S111 > Rmax
    S111r = Rmax;
elseif S111 < Rmin
    S111r = Rmin;
end
end

%% realizability_S2
function [S110r,S101r,S011r,S2r] = realizability_S2(S110,S101,S011)
% realizability_S2 checks and corrects realizablity of 2nd-order moments
%
S2 = 1 + 2*S110*S101*S011 - (S110^2+S101^2+S011^2);
xr = 1;
if S2 < 0
    Y = @(x) 1 + 2*S110*S101*S011*x.^3 - (S110^2+S101^2+S011^2)*x.^2;
    xr = fzero(Y,[0 1]);
end
xr = 0.9999*xr;
S110r = xr*S110;
S101r = xr*S101;
S011r = xr*S011;
S2r = 1 + 2*S110r*S101r*S011r - (S110r^2+S101r^2+S011r^2);
if S2r < 0
    warning('S2 < 0 after correction in realizability_S2')
    disp([S2,S2r])
end
end

%% realizability_S210
function [S210r,S201r] = realizability_S210(S110,S101,S011,S300,S210,S201,H200,beta)
% realizability_S210 checks and corrects realizablity of S210, S201
%
xr = 1;
X = [S210-S110*S300; S201-S101*S300];
D1 = [1-S101^2, S101*S110-S011; S101*S110-S011, 1-S110^2];
U = sqrtm(D1);
V = U*X;
L = max([0 V'*V]);
dD1 = max([0 det(D1)]);
R = H200*dD1;
if R <= 0 || X'*X < 1000*eps
    xr = 0;
elseif L > R
    xr = sqrt(R/L);
end
Vr = beta*xr*V;
Xr = U\Vr;
S210r = Xr(1) + S110*S300;
S201r = Xr(2) + S101*S300;
end

%% realizability_S211
function S211r = realizability_S211(e11,e22,e33,e12,e13,d23,S211,beta)
% realizability_S211
 
S211r = S211;
b211 = e12*e13;
G211 = max([0 (e11*e22-e13^2)*(e11*e33-e12^2)]);
sG211 = beta*sqrt(G211);
s211min = d23+(b211-sG211)/e11;
s211max = d23+(b211+sG211)/e11;
%
if S211 <= s211min
    S211r = s211min;
elseif S211 >= s211max
    S211r = s211max;
end
%
end

%% realizability_S310
function [S310r,S220r] = realizability_S310(S110,S101,S011,S300,S210,S201,S120,S111,S310,S220,H200,beta)

% realizability_S310 checks and corrects realizablity of S310 and S220
 
S310r = S310;
S220r = S220;
 
b310 = S111*((1 - S110^2)*S201 + (S101*S110 - S011)*S210 + (S011*S110 - S101)*S300) ...
     + S120*((1 - S101^2)*S210 + (S011*S101 - S110)*S300 + (S101*S110 - S011)*S201) ...
     + S210*((1 - S011^2)*S300 + (S011*S101 - S110)*S210 + (S011*S110 - S101)*S201);
 
D1 = [1-S101^2, S101*S110-S011;S101*S110-S011, 1-S110^2];
dD1 = det(D1);
V1 = [S210-S110*S300; S201-S101*S300];
L1 = V1'*D1*V1/dD1;
 
G310b = H200 - L1;
if G310b < 0
    G310b = 0;
end
 
D2 = [1-S011^2, S011*S101-S110, S011*S110-S101;...
      S011*S101-S110, 1-S101^2, S101*S110-S011;...
      S011*S110-S101, S101*S110-S011, 1-S110^2];
V2 = [S210; S120; S111];
L2 = V2'*D2*V2/dD1;
 
G310a = S220 - S110^2 - L2 + 1000*eps;
if G310a < 0 || L2 < 0
    G310a = 0;
end

G310 = G310a*G310b;
sG310 = beta*sqrt(G310);

s310min = S110 + b310/dD1 - sG310;
s310max = S110 + b310/dD1 + sG310;

if S310 <= s310min
    S310r = s310min;
elseif S310 >= s310max
    S310r = s310max;
end
%
end

%% realizability_S310_220
function [S220r] = realizability_S310_220(S110,S101,S011,S210,S120,S111,S220)
% realizability_S310_220 checks and corrects realizablity of S220 for S310

S220r = S220;
D1 = [1-S101^2, S101*S110-S011; S101*S110-S011, 1-S110^2];
dD1 = det(D1);
 
D2 = [1-S011^2, S011*S101-S110, S011*S110-S101;...
      S011*S101-S110, 1-S101^2, S101*S110-S011;...
      S011*S110-S101, S101*S110-S011, 1-S110^2];
V2 = [S210; S120; S111];
L2 = V2'*D2*V2/dD1;
S220 = S110^2 + L2 + 1000*eps;
 
if S220 > S220r
    S220r = S220;
end
 
end

%% realizablity_S220
function S220r = realizablity_S220(S110,S220,A220)
% realizablity_S220 checks maximum bounds and corrects S220

S220r = S220;
s220min = max([S110^2 1-A220]);
s220max = 1+A220;
if S220 < s220min
    S220r = s220min;
elseif S220 > s220max
    S220r = s220max;
end
end

