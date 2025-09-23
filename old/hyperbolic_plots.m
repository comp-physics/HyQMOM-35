% contour plots of eigenvalues of J6 to check for hyperbolicity
figure(12)
cmin = 0; cmax = 100;

lam6x = zeros(Np,Np,6);
lam6y = zeros(Np,Np,6);
for i = 1:Np
    for j = 1:Np
        M1 = zeros(Nmom,1);
        for k=1:Nmom
            M1(k,1) = M(i,j,k);
        end
        %
        m00 = M1(1);
        m10 = M1(2);
        m20 = M1(3);
        m30 = M1(4);
        m40 = M1(5);
        m01 = M1(6);
        m11 = M1(7);
        m21 = M1(8);
        m31 = M1(9);
        m02 = M1(10);
        m12 = M1(11);
        m22 = M1(12);
        m03 = M1(13);
        m13 = M1(14);
        m04 = M1(15);
        J6x = jacobian6(m00,m01,m02,m03,m04,m10,m11,m12,m13,m20,m21,m22,m30,m31,m40);
        lam6x(i,j,:) = eig(J6x);
        lam6x(i,j,:) = sort(real(lam6x(i,j,:)));
        % J6y = jacobian6(m00,m10,m20,m30,m40,m01,m11,m21,m31,m02,m12,m22,m03,m13,m04);
        % lam6y(i,j,:) = eig(J6y);
        % lam6y(i,j,:) = sort(real(lam6y(i,j,:)));
    end
end

subplot(3,3,1)
[X,Y] = meshgrid(xm,ym);
Z = lam6x(:,:,1);
contourf(X,Y,Z',50)
axis square
title('\lambda_1')

subplot(3,3,2)
[X,Y] = meshgrid(xm,ym);
Z = lam6x(:,:,2);
contourf(X,Y,Z',50)
axis square
title('\lambda_2')

subplot(3,3,3)
[X,Y] = meshgrid(xm,ym);
Z = lam6x(:,:,3);
contourf(X,Y,Z',50)
axis square
title('\lambda_3')

subplot(3,3,4)
[X,Y] = meshgrid(xm,ym);
Z = lam6x(:,:,4);
contourf(X,Y,Z',50)
axis square
title('\lambda_4')

subplot(3,3,5)
[X,Y] = meshgrid(xm,ym);
Z = lam6x(:,:,5);
contourf(X,Y,Z',50)
axis square
title('\lambda_5')

subplot(3,3,6)
[X,Y] = meshgrid(xm,ym);
Z = lam6x(:,:,6);
contourf(X,Y,Z',50)
axis square
title('\lambda_6')

subplot(3,3,7)
[X,Y] = meshgrid(xm,ym);
Z = lam6x(:,:,2) - lam6x(:,:,1);
contourf(X,Y,Z',50)
axis square
title('\lambda_2-\lambda_1')

subplot(3,3,8)
[X,Y] = meshgrid(xm,ym);
Z = lam6x(:,:,4) - lam6x(:,:,3);
contourf(X,Y,Z',50)
axis square
title('\lambda_4-\lambda_3')

subplot(3,3,9)
[X,Y] = meshgrid(xm,ym);
Z = lam6x(:,:,6) - lam6x(:,:,5);
contourf(X,Y,Z',50)
axis square
title('\lambda_6-\lambda_5')