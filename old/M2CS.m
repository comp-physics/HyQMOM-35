function [C,S] = M2CS(M)
% M(1:Nmom+1): moments of order 0 to Nmom-1
% 
% normalized central moments
    Nmom = size(M,2);
    C = zeros(1,Nmom-1);
    S = zeros(1,Nmom-1);
    Mn = M(2:end)/M(1);
    rho = M(1);
    u = Mn(1);
    C(2) = Mn(2) - u^2;
    p = C(2);
    S(2) = 1;
    for k = 3:Nmom-1
        C(k) = (-1)^k*u^k;
        for j = 1:k
            C(k) = C(k) + nchoosek(k,j)*(-u)^(k-j)*Mn(j);
        end
        S(k) = C(k)/p^(k/2);
    end
end