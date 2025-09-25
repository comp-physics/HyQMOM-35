function varargout = realizability_engine(type, varargin)
% realizability_engine - Unified realizability checks and corrections
%
% Usage:
%   [S110r,S101r,S011r,S2r] = realizability_engine('S2', S110, S101, S011)
%   [S111r] = realizability_engine('S111', S110, S101, S011, S210, S201, S120, S021, S102, S012, S111)
%   [S210r,S201r] = realizability_engine('S210', S110, S101, S011, S300, S210, S201, H200, beta)
%   [S211r] = realizability_engine('S211', e11, e22, e33, e12, e13, d23, S211, beta)
%   [S220r] = realizability_engine('S220', S110, S220, A220)
%   [S310r,S220r] = realizability_engine('S310', S110, S101, S011, S300, S210, S201, S120, S111, S310, S220, H200, beta)
%   [S220r] = realizability_engine('S310_220', S110, S210, S120, S111, A220, G310a, S220)

switch upper(type)
    case 'S2'
        % realizability_S2 checks and corrects realizability of 2nd-order moments
        S110 = varargin{1};
        S101 = varargin{2};
        S011 = varargin{3};
        
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
        varargout = {S110r, S101r, S011r, S2r};
        
    case 'S111'
        % realizability_S111 checks and corrects realizability of S111
        S110 = varargin{1};
        S101 = varargin{2};
        S011 = varargin{3};
        S210 = varargin{4};
        S201 = varargin{5};
        S120 = varargin{6};
        S021 = varargin{7};
        S102 = varargin{8};
        S012 = varargin{9};
        S111 = varargin{10};
        
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
        varargout = {S111r};
        
    case 'S210'
        % realizability_S210 checks and corrects realizability of S210, S201
        S110 = varargin{1};
        S101 = varargin{2};
        S011 = varargin{3};
        S300 = varargin{4};
        S210 = varargin{5};
        S201 = varargin{6};
        H200 = varargin{7};
        beta = varargin{8};
        
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
        varargout = {S210r, S201r};
        
    case 'S211'
        % realizability_S211 checks and corrects realizability of S211
        e11 = varargin{1};
        e22 = varargin{2};
        e33 = varargin{3};
        e12 = varargin{4};
        e13 = varargin{5};
        d23 = varargin{6};
        S211 = varargin{7};
        beta = varargin{8};
        
        S211r = S211;
        b211 = e12*e13;
        G211 = max([0 (e11*e22-e13^2)*(e11*e33-e12^2)]);
        sG211 = beta*sqrt(G211);
        s211min = d23+(b211-sG211)/e11;
        s211max = d23+(b211+sG211)/e11;
        
        if S211 <= s211min
            S211r = s211min;
        elseif S211 >= s211max
            S211r = s211max;
        end
        
        varargout = {S211r};
        
    case 'S220'
        % realizability_S220 checks realizability of S220
        S110 = varargin{1};
        S220 = varargin{2};
        A220 = varargin{3};
        
        s220max = A220 - S110^2;
        S220r = min(S220, s220max);
        varargout = {S220r};
        
    case 'S310'
        % realizability_S310 checks and corrects realizability of S310 and S220
        S110 = varargin{1};
        S101 = varargin{2};
        S011 = varargin{3};
        S300 = varargin{4};
        S210 = varargin{5};
        S201 = varargin{6};
        S120 = varargin{7};
        S111 = varargin{8};
        S310 = varargin{9};
        S220 = varargin{10};
        H200 = varargin{11};
        beta = varargin{12};
        
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
            warning('G310b < 0')
            G310b = 0;
        end
        
        D2 = [1-S011^2, S011*S101-S110, S011*S110-S101;...
              S011*S101-S110, 1-S101^2, S101*S110-S011;...
              S011*S110-S101, S101*S110-S011, 1-S110^2];
        V2 = [S210; S120; S111];
        L2 = V2'*D2*V2/dD1;
        
        G310a = S220 - S110^2 - L2 + 1000*eps;
        if G310a < 0 || L2 < 0
            warning('G310a < 0 || L2 < 0')
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
        
        varargout = {S310r, S220r};
        
    case 'S310_220'
        % realizability_S310_220 returns realizability corrected S220
        S110 = varargin{1};
        S210 = varargin{2};
        S120 = varargin{3};
        S111 = varargin{4};
        A220 = varargin{5};
        G310a = varargin{6};
        S220 = varargin{7};
        
        D2 = [S210; S120; S111];
        L2 = D2'*D2;
        if G310a >= eps
            S220r = S110^2 + L2 + 1000*eps;
        else
            s22min = S110^2 + L2 + 1000*eps;
            S220r = max(S220,s22min);
            S220r = min(S220r,A220-S110^2);
        end
        varargout = {S220r};
        
    otherwise
        error('realizability_engine: Unknown type %s', type);
end

end
