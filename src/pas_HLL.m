function Mp = pas_HLL(M,F,dt,dx,vpmin,vpmax)
    Np = size(M,1);
    Nmom = size(M,2);
    
    % Vectorized edge speed computation
    lleft = min(vpmin(1:end-1), vpmin(2:end));
    lright = max(vpmax(1:end-1), vpmax(2:end));
    denom = (lleft - lright);
    safe = abs(denom) > 1e-10;
    
    % Vectorized star state computation
    Wstar = zeros(Np,Nmom);
    j_safe = find(safe);
    for k = 1:length(j_safe)
        j = j_safe(k);
        Wstar(j,:) = (lleft(j)*M(j,:) - lright(j)*M(j+1,:))/denom(j) ...
                   - (F(j,:) - F(j+1,:))/denom(j);
    end
    
    % Boundary conditions for Wstar
    Wstar(1,:) = Wstar(2,:);
    Wstar(Np,:) = Wstar(Np-1,:);
    
    % Boundary conditions for F
    F(1,:) = F(2,:);
    F(Np,:) = F(Np-1,:);
    
    % Vectorized flux computation (inline flux_HLL)
    Flux = zeros(Np-1, Nmom);
    for j = 1:Np-1
        Flux(j,:) = 0.5*(F(j,:) + F(j+1,:)) ...
                  - 0.5*(abs(lleft(j))*(Wstar(j,:) - M(j,:)) ...
                       - abs(lright(j))*(Wstar(j,:) - M(j+1,:)));
    end
    
    % Time stepping
    Mp = M;
    Mp(2:Np-1,:) = M(2:Np-1,:) - (dt/dx)*(Flux(2:Np-1,:) - Flux(1:Np-2,:));
    
    % Boundary conditions
    Mp(1,:) = Mp(2,:);
    Mp(Np,:) = Mp(Np-1,:);
    
end