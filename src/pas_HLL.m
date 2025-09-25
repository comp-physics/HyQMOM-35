function Mp = pas_HLL(M,F,dt,dx,vpmin,vpmax)
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
    Wstar(1,:)  = Wstar(2,:);
    Wstar(Np,:) = Wstar(Np-1,:);
    F(1,:)  = F(2,:);
    F(Np,:) = F(Np-1,:);

    Flux = flux_HLL(Wstar,M,lleft,lright,F,Np);

    Mp = M;
    Mp(2:Np-1,:) = M(2:Np-1,:) - dt/dx*(Flux(2:Np-1,:)-Flux(1:Np-2,:));

    %BC
    Mp(1,:) = Mp(2,:);
    Mp(Np,:) = Mp(Np-1,:);

    % figure(20)
    % plot(lleft,'o')
    % hold on
    % plot(lright,'p')

end