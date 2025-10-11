% contour plots of eigenvalues of J6 to check for hyperbolicity
figure(12)
cmin = 0; cmax = 100;

lam6x = zeros(Np,Np,6);
lam6y = zeros(Np,Np,6);
lam6z = zeros(Np,Np,6);
for i = 1:Np
    for j = 1:Np
        M1 = zeros(Nmom,1);
        for k=1:Nmom
            M1(k,1) = M(i,j,k);
        end
        %
        m000 = M1(1);
        m100 = M1(2);
        m200 = M1(3);
        m300 = M1(4);
        m400 = M1(5);
        m010 = M1(6);
        m110 = M1(7);
        m210 = M1(8);
        m310 = M1(9);
        m020 = M1(10);
        m120 = M1(11);
        m220 = M1(12);
        m030 = M1(13);
        m130 = M1(14);
        m040 = M1(15);
        m001 = M1(16);
        m101 = M1(17);
        m201 = M1(18);
        m301 = M1(19);
        m002 = M1(20);
        m102 = M1(21);
        m202 = M1(22);
        m003 = M1(23);
        m103 = M1(24);
        m004 = M1(25);
        m011 = M1(26);
        m021 = M1(29);
        m031 = M1(31);
        m012 = M1(32);
        m013 = M1(34);
        m022 = M1(35);
        % UV  moments
        J6 = jacobian6(m000,m010,m020,m030,m040,m100,m110,m120,m130,m200,m210,m220,m300,m310,m400);
        lam6x(i,j,:) = eig(J6);
        lam6x(i,j,:) = sort(real(lam6x(i,j,:)));
        % UW  moments
        J6 = jacobian6(m000,m001,m002,m003,m004,m100,m101,m102,m103,m200,m201,m202,m300,m301,m400);
        lam6y(i,j,:) = eig(J6);
        lam6y(i,j,:) = sort(real(lam6y(i,j,:)));
        % VW  moments
        J6 = jacobian6(m000,m001,m002,m003,m004,m010,m011,m012,m013,m020,m021,m022,m030,m031,m040);
        lam6z(i,j,:) = eig(J6);
        lam6z(i,j,:) = sort(real(lam6z(i,j,:)));
    end
end

subplot(3,3,1)
[X,Y] = meshgrid(xm,ym);
Z = lam6x(:,:,1);
contourf(X,Y,Z',50)
axis square
title('\lambda_1')
colorbar

subplot(3,3,2)
[X,Y] = meshgrid(xm,ym);
Z = lam6x(:,:,2);
contourf(X,Y,Z',50)
axis square
title('\lambda_2')
colorbar

subplot(3,3,3)
[X,Y] = meshgrid(xm,ym);
Z = lam6x(:,:,3);
contourf(X,Y,Z',50)
axis square
title('\lambda_3')
colorbar

subplot(3,3,4)
[X,Y] = meshgrid(xm,ym);
Z = lam6x(:,:,4);
contourf(X,Y,Z',50)
axis square
title('\lambda_4')
colorbar

subplot(3,3,5)
[X,Y] = meshgrid(xm,ym);
Z = lam6x(:,:,5);
contourf(X,Y,Z',50)
axis square
title('\lambda_5')
colorbar

subplot(3,3,6)
[X,Y] = meshgrid(xm,ym);
Z = lam6x(:,:,6);
contourf(X,Y,Z',50)
axis square
title('\lambda_6')
colorbar

subplot(3,3,7)
[X,Y] = meshgrid(xm,ym);
Z = lam6x(:,:,2) - lam6x(:,:,1);
contourf(X,Y,Z',50)
axis square
title('\lambda_2-\lambda_1')
colorbar

subplot(3,3,8)
[X,Y] = meshgrid(xm,ym);
Z = lam6x(:,:,4) - lam6x(:,:,3);
contourf(X,Y,Z',50)
axis square
title('\lambda_4-\lambda_3')
colorbar

subplot(3,3,9)
[X,Y] = meshgrid(xm,ym);
Z = lam6x(:,:,6) - lam6x(:,:,5);
contourf(X,Y,Z',50)
axis square
title('\lambda_6-\lambda_5')
colorbar

colormap sky

figure(13)

subplot(3,3,1)
[X,Y] = meshgrid(xm,ym);
Z = lam6y(:,:,1);
contourf(X,Y,Z',50)
axis square
title('\lambda_1')
colorbar

subplot(3,3,2)
[X,Y] = meshgrid(xm,ym);
Z = lam6y(:,:,2);
contourf(X,Y,Z',50)
axis square
title('\lambda_2')
colorbar

subplot(3,3,3)
[X,Y] = meshgrid(xm,ym);
Z = lam6y(:,:,3);
contourf(X,Y,Z',50)
axis square
title('\lambda_3')
colorbar

subplot(3,3,4)
[X,Y] = meshgrid(xm,ym);
Z = lam6y(:,:,4);
contourf(X,Y,Z',50)
axis square
title('\lambda_4')
colorbar

subplot(3,3,5)
[X,Y] = meshgrid(xm,ym);
Z = lam6y(:,:,5);
contourf(X,Y,Z',50)
axis square
title('\lambda_5')
colorbar

subplot(3,3,6)
[X,Y] = meshgrid(xm,ym);
Z = lam6y(:,:,6);
contourf(X,Y,Z',50)
axis square
title('\lambda_6')
colorbar

subplot(3,3,7)
[X,Y] = meshgrid(xm,ym);
Z = lam6y(:,:,2) - lam6y(:,:,1);
contourf(X,Y,Z',50)
axis square
title('\lambda_2-\lambda_1')
colorbar

subplot(3,3,8)
[X,Y] = meshgrid(xm,ym);
Z = lam6y(:,:,4) - lam6y(:,:,3);
contourf(X,Y,Z',50)
axis square
title('\lambda_4-\lambda_3')
colorbar

subplot(3,3,9)
[X,Y] = meshgrid(xm,ym);
Z = lam6y(:,:,6) - lam6y(:,:,5);
contourf(X,Y,Z',50)
axis square
title('\lambda_6-\lambda_5')
colorbar

colormap sky

figure(14)

subplot(3,3,1)
[X,Y] = meshgrid(xm,ym);
Z = lam6z(:,:,1);
contourf(X,Y,Z',50)
axis square
title('\lambda_1')
colorbar

subplot(3,3,2)
[X,Y] = meshgrid(xm,ym);
Z = lam6z(:,:,2);
contourf(X,Y,Z',50)
axis square
title('\lambda_2')
colorbar

subplot(3,3,3)
[X,Y] = meshgrid(xm,ym);
Z = lam6z(:,:,3);
contourf(X,Y,Z',50)
axis square
title('\lambda_3')
colorbar

subplot(3,3,4)
[X,Y] = meshgrid(xm,ym);
Z = lam6z(:,:,4);
contourf(X,Y,Z',50)
axis square
title('\lambda_4')
colorbar

subplot(3,3,5)
[X,Y] = meshgrid(xm,ym);
Z = lam6z(:,:,5);
contourf(X,Y,Z',50)
axis square
title('\lambda_5')
colorbar

subplot(3,3,6)
[X,Y] = meshgrid(xm,ym);
Z = lam6z(:,:,6);
contourf(X,Y,Z',50)
axis square
title('\lambda_6')
colorbar

subplot(3,3,7)
[X,Y] = meshgrid(xm,ym);
Z = lam6z(:,:,2) - lam6z(:,:,1);
contourf(X,Y,Z',50)
axis square
title('\lambda_2-\lambda_1')
colorbar

subplot(3,3,8)
[X,Y] = meshgrid(xm,ym);
Z = lam6z(:,:,4) - lam6z(:,:,3);
contourf(X,Y,Z',50)
axis square
title('\lambda_4-\lambda_3')
colorbar

subplot(3,3,9)
[X,Y] = meshgrid(xm,ym);
Z = lam6z(:,:,6) - lam6z(:,:,5);
contourf(X,Y,Z',50)
axis square
title('\lambda_6-\lambda_5')
colorbar

colormap sky