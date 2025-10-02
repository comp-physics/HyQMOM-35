function Mp = pas_HLL(M,F,dt,dx,vpmin,vpmax,apply_bc_left,apply_bc_right)
% pas_HLL - HLL flux update
%
% Optional inputs:
%   apply_bc_left/right - (optional) Apply BCs at boundaries (default: true for both)
%                         Set to false for processor boundaries in MPI

    if nargin < 7
        apply_bc_left = true;
    end
    if nargin < 8
        apply_bc_right = true;
    end

    Np = size(M,1);
    Nmom = size(M,2);
    Wstar = zeros(Np,Nmom);
    lleft = zeros(Np,1);
    lright = zeros(Np,1);
    
    % Determine loop range based on boundary conditions
    % If left processor boundary (apply_bc_left=false), extend loop to START at 1
    j_start = 2;
    if ~apply_bc_left
        j_start = 1;  % Start at 1 to compute Wstar(1) using left neighbor M(1)
    end
    
    % For right boundary, we always end at Np-1 because computing at Np would need M(Np+1)
    j_end = Np-1;
    
    % Compute Wstar using stencil (needs M(j) and M(j+1))
    for j = j_start:j_end
        lleft(j) = min([vpmin(j),vpmin(j+1)]);
        lright(j) = max([vpmax(j),vpmax(j+1)]);
        %Wstar
        if abs(lleft(j)-lright(j))>1.d-10
            Wstar(j,:) = (lleft(j)*M(j,:) - lright(j)*M(j+1,:))/(lleft(j) - lright(j))...
                - (F(j,:)-F(j+1,:))/(lleft(j)-lright(j));
        else
            Wstar(j,:)=0.;
        end
    end
    
    % Apply boundary conditions ONLY at physical boundaries
    if apply_bc_left
        Wstar(1,:) = Wstar(2,:);
    end
    if apply_bc_right
        Wstar(Np,:) = Wstar(Np-1,:);
    end
    
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