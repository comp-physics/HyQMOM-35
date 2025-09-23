% contour plots of various quantities 

cmin = 0; cmax = 100;
for i = 1:Np
    for j = 1:Np
        s30 = S(i,j,4);
        s40 = S(i,j,5);
        s11 = S(i,j,7);
        s21 = S(i,j,8);
        s31 = S(i,j,9);
        s12 = S(i,j,11);
        s22 = S(i,j,12);
        s03 = S(i,j,13);
        s13 = S(i,j,14);
        s04 = S(i,j,15);
        [D2,~] = delta2star(s03,s04,s11,s12,s13,s21,s22,s30,s31,s40);
        lamdaD2 = eig(D2);
        lamdaD2 = sort(real(lamdaD2));
        LAM1D2(i,j) = min(cmax,max(cmin,lamdaD2(1)));
        LAM2D2(i,j) = min(cmax,max(cmin,lamdaD2(2)));
        LAM3D2(i,j) = min(cmax,max(cmin,lamdaD2(3)));
    end
end

subplot(3,3,1)
[X,Y] = meshgrid(xm,ym);
Z = LAM1D2;
contourf(X,Y,Z',50)
axis square
title('\lambda_1^*')

subplot(3,3,2)
[X,Y] = meshgrid(xm,ym);
Z = LAM2D2;
contourf(X,Y,Z',50)
axis square
title('\lambda_2^*')

subplot(3,3,3)
[X,Y] = meshgrid(xm,ym);
Z = LAM3D2;
contourf(X,Y,Z',50)
axis square
title('\lambda_3^*')

for i = 1:Np
    for j = 1:Np
        s30 = S(i,j,4);
        s40 = S(i,j,5);
        s11 = S(i,j,7);
        s21 = S(i,j,8);
        s31 = S(i,j,9);
        s12 = S(i,j,11);
        s22 = S(i,j,12);
        s03 = S(i,j,13);
        s13 = S(i,j,14);
        s04 = S(i,j,15);
        [s21,s12,s31,s22,s13] = realizable_2D(s30,s40,s11,s21,s31,s12,s22,s03,s13,s04);
        Lambda = delta2starchol(s03,s04,s11,s12,s13,s21,s22,s30,s31,s40);
        lamdaD2 = Lambda;
        LAM1D2(i,j) = min(cmax,max(cmin,lamdaD2(1)));
        LAM2D2(i,j) = min(cmax,max(cmin,lamdaD2(2)));
        LAM3D2(i,j) = min(cmax,max(cmin,lamdaD2(3)));
        %
        Del1= max(eps,1 - s11^2);
        H20 = max(eps,s40 - s30^2 - 1);
        H02 = max(eps,s04 - s03^2 - 1);
        R = rootsR(Del1,H02,H20,s03,s11,s12,s13,s21,s30,s31);
        Rr = sort(real(R));
        Root1(i,j) = Rr(1);
        Root2(i,j) = Rr(2);
        Root3(i,j) = Rr(3);
    end
end

subplot(3,3,4)
[X,Y] = meshgrid(xm,ym);
Z = LAM1D2;
contourf(X,Y,Z',50)
axis square
title('\gamma_1^*')

subplot(3,3,5)
[X,Y] = meshgrid(xm,ym);
Z = LAM2D2;
contourf(X,Y,Z',50)
axis square
title('\gamma_2^*')

subplot(3,3,6)
[X,Y] = meshgrid(xm,ym);
Z = LAM3D2;
contourf(X,Y,Z',50)
axis square
title('\gamma_3^*')

subplot(3,3,7)
[X,Y] = meshgrid(xm,ym);
Z = Root1;
contourf(X,Y,Z',50)
axis square
title('R_{1}')

subplot(3,3,8)
[X,Y] = meshgrid(xm,ym);
Z = Root2;
contourf(X,Y,Z',50)
axis square
title('R_{2}')

subplot(3,3,9)
[X,Y] = meshgrid(xm,ym);
Z = Root3;
contourf(X,Y,Z',50)
axis square
title('R_{3}')