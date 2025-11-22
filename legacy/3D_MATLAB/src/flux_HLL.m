function Flux = flux_HLL(Wstar,W,l1,l2,F,N)
for j = 1:N-1
    Flux(j,:) = 1/2*(F(j,:) + F(j+1,:)) - 1/2*( abs(l1(j))*(Wstar(j,:) - W(j,:)) - abs(l2(j))*(Wstar(j,:) - W(j+1,:)) );
end

end