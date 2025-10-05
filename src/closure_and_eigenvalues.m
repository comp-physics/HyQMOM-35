function [Mp,vpmin,vpmax] = closure_and_eigenvalues(mom)
    % mom(1:2N+1): moments of order 0 to 2N
    % Mp: moment of order 2N+1

    N = (size(mom,2)-1)/2;

    % recurrence coefficients
    sig = zeros(N+2,2*N+3);
    a = zeros(N+1,1) ;
    b = zeros(N+2,1) ;
    for i = 2:(2*N+2)
        sig(2,i) = mom(i-1) ;
    end
    a(1) = mom(2)/mom(1) ;
    for k = 3:(N+1)
        for l = k:(2*N-k+4)
            sig(k,l) = sig(k-1,l+1)-a(k-2)*sig(k-1,l)-b(k-2)*sig(k-2,l) ;
        end
        a(k-1) = sig(k,k+1)/sig(k,k)-sig(k-1,k)/sig(k-1,k-1) ;
        b(k-1) = sig(k,k)/sig(k-1,k-1) ;
    end
    k = N+2;
    sig(k,k) = sig(k-1,k+1)-a(k-2)*sig(k-1,k)-b(k-2)*sig(k-2,k) ;
    b(k-1) = sig(k,k)/sig(k-1,k-1) ;

    % closure
    a(N+1) = sum(a(1:N))/N;

    % moment of order 2N+1
    sig(N+2,N+3) = sig(N+2,N+2)*(a(N+1)+sig(N+1,N+2)/sig(N+1,N+1));
    for k = N+2:-1:3
        l = 2*N-k+5;
        sig(k-1,l+1) = sig(k,l) +a(k-2)*sig(k-1,l)+b(k-2)*sig(k-2,l) ;
    end
    Mp = sig(2,2*N+3);
    if any(b<0)
        disp('pb realizabilite')
        b
        mom
    end

    % Computation of the maximal and minimal values of the eigenvalues
    
    % Setup Jacobi matrix to find roots of R_{n+1}
    b(N+1) = b(N+1)*(2*N+1)/N;
    z = zeros(N+1,N+1) ;
    for i = 1:N
        z(i,i) = a(i) ;
        z(i,i+1) = sqrt(b(i+1)) ;
        z(i+1,i) = z(i,i+1) ;
    end
    z(N+1,N+1) = a(N+1) ;

    % Compute weights and abscissas
    vp = eig(z) ;
    vpmin = min(vp);
    vpmax = max(vp);
end