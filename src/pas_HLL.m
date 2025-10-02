function Mp = pas_HLL(M,F,dt,dx,vpmin,vpmax,apply_bc_left,apply_bc_right)
% pas_HLL - HLL flux update with optional boundary condition control
%
% Optional inputs:
%   apply_bc_left  - (optional) true to apply BC at left boundary (default: true)
%   apply_bc_right - (optional) true to apply BC at right boundary (default: true)
%
% For MPI: set to false for processor boundaries, true for physical boundaries

    if nargin < 7
        apply_bc_left = true;
    end
    if nargin < 8
        apply_bc_right = true;
    end

    Np = size(M,1);
    Nmom = size(M,2);
    Wstar = zeros(Np,Nmom);
    for j = 2:Np-1
        lleft(j) =min([vpmin(j),vpmin(j+1)]);
        lright(j)=max([vpmax(j),vpmax(j+1)]);
        %Wstar
        if abs(lleft(j)-lright(j))>1.d-10
            Wstar(j,:) = (lleft(j)*M(j,:) - lright(j)*M(j+1,:))/(lleft(j) - lright(j))...
                - (F(j,:)-F(j+1,:))/(lleft(j)-lright(j));
        else
            Wstar(j,:)=0.;
        end
    end
    
    % Apply boundary conditions based on boundary type
    Wstar(1,:)  = Wstar(2,:);
    Wstar(Np,:) = Wstar(Np-1,:);
    
    % Only apply flux BCs at physical boundaries
    if apply_bc_left
        F(1,:) = F(2,:);
    end
    if apply_bc_right
        F(Np,:) = F(Np-1,:);
    end

    Flux = flux_HLL(Wstar,M,lleft,lright,F,Np);

    Mp = M;
    Mp(2:Np-1,:) = M(2:Np-1,:) - dt/dx*(Flux(2:Np-1,:)-Flux(1:Np-2,:));

    % Apply solution BCs only at physical boundaries
    if apply_bc_left
        Mp(1,:) = Mp(2,:);
    end
    if apply_bc_right
        Mp(Np,:) = Mp(Np-1,:);
    end

    % figure(20)
    % plot(lleft,'o')
    % hold on
    % plot(lright,'p')

end