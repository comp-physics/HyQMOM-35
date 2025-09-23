function M = S2M(S,rho,u,p)
% 
    Nmom = size(S,2)+1;
    M = zeros(1,Nmom);
    C = zeros(1,Nmom-1);
    Mn = zeros(1,Nmom-1);
    M(1) = rho;
    Mn(1) = u;
    Mn(2) = u^2+p;
    C(2) = p;
    for k = 3:Nmom-1
        C(k) = S(k)*p^(k/2);
        Mn(k) = u^k;
        for j = 2:k
            Mn(k) = Mn(k) + nchoosek(k,j)*u^(k-j)*C(j);
        end
    end
    M(2:end) = rho*Mn;
end