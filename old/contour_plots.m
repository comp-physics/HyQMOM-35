% contour plots of various quantities 

subplot(3,3,1)
[X,Y] = meshgrid(xm,ym);
Z = M(:,:,1);
contourf(X,Y,Z',50)
axis square
title('Density')

subplot(3,3,2)
[X,Y] = meshgrid(xm,ym);
Z = M(:,:,2)./M(:,:,1);
contourf(X,Y,Z',50)
axis square
title('U velocity')

subplot(3,3,3)
[X,Y] = meshgrid(xm,ym);
Z = M(:,:,6)./M(:,:,1);
contourf(X,Y,Z',50)
axis square
title('V velocity')

subplot(3,3,4)
[X,Y] = meshgrid(xm,ym);
Z = C(:,:,3);
contourf(X,Y,Z',50)
axis square
title('C_{20}')

subplot(3,3,5)
[X,Y] = meshgrid(xm,ym);
Z = C(:,:,10);
contourf(X,Y,Z',50)
axis square
title('C_{02}')

subplot(3,3,6)
[X,Y] = meshgrid(xm,ym);
Z = 1 - S(:,:,7).^2;
%Z = max(0,Z);
contourf(X,Y,Z',50)
axis square
title('\Delta_1')

subplot(3,3,7)
[X,Y] = meshgrid(xm,ym);
Z = max(0,S(:,:,5) - S(:,:,4).^2 - 1);
contourf(X,Y,Z',50)
axis square
title('H_{20}')

subplot(3,3,8)
[X,Y] = meshgrid(xm,ym);
Z = max(0,S(:,:,15) - S(:,:,13).^2 - 1);
contourf(X,Y,Z',50)
axis square
title('H_{02}')

subplot(3,3,9)
[X,Y] = meshgrid(xm,ym);
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
        dD2s = det(D2);
        if dD2s >= -10
            dDel2(i,j) = dD2s;
        else
            dDel2(i,j) = -10;
        end
    end
end
Z = dDel2;
contourf(X,Y,Z',50)
axis square
title('\Delta_2^*')