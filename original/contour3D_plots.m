% contour plots of various quantities 
% figure(9)
subplot(3,4,1)
[X,Y] = meshgrid(xm,ym);
Z = M(:,:,1);
contourf(X,Y,Z',50,'k')
axis square
title('Density')
colorbar

subplot(3,4,2)
[X,Y] = meshgrid(xm,ym);
Z = M(:,:,2)./M(:,:,1);
contourf(X,Y,Z',50,'k')
axis square
title('U velocity')
colorbar

subplot(3,4,3)
[X,Y] = meshgrid(xm,ym);
Z = M(:,:,6)./M(:,:,1);
contourf(X,Y,Z',50)
axis square
title('V velocity')
colorbar

subplot(3,4,4)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,7);
contourf(X,Y,Z',50)
axis square
title('S_{110}')
colorbar

subplot(3,4,5)
[X,Y] = meshgrid(xm,ym);
Z = C(:,:,3);
contourf(X,Y,Z',50)
axis square
title('C_{200}')
colorbar

subplot(3,4,6)
[X,Y] = meshgrid(xm,ym);
Z = C(:,:,10);
contourf(X,Y,Z',50)
axis square
title('C_{020}')
colorbar

subplot(3,4,7)
[X,Y] = meshgrid(xm,ym);
Z = C(:,:,20);
contourf(X,Y,Z',50)
axis square
title('C_{002}')
colorbar

subplot(3,4,8)
[X,Y] = meshgrid(xm,ym);
Z = 1 + 2*S(:,:,7).*S(:,:,17).*S(:,:,26) - S(:,:,7).^2 - S(:,:,17).^2 - S(:,:,26).^2;
%Z = max(0,Z);
contourf(X,Y,Z',50)
axis square
title('|\Delta_1|')
colorbar

subplot(3,4,9)
[X,Y] = meshgrid(xm,ym);
Z = max(0,S(:,:,5) - S(:,:,4).^2 - 1);
contourf(X,Y,Z',50)
axis square
title('H_{200}')
colorbar

subplot(3,4,10)
[X,Y] = meshgrid(xm,ym);
Z = max(0,S(:,:,15) - S(:,:,13).^2 - 1);
contourf(X,Y,Z',50)
axis square
title('H_{020}')
colorbar

subplot(3,4,11)
[X,Y] = meshgrid(xm,ym);
Z = max(0,S(:,:,25) - S(:,:,23).^2 - 1);
contourf(X,Y,Z',50)
axis square
title('H_{002}')
colorbar

subplot(3,4,12)
[X,Y] = meshgrid(xm,ym);
for i = 1:Np
    for j = 1:Np
        S300 = S(i,j,4);
        S400 = S(i,j,5);
        S110 = S(i,j,7);
        S210 = S(i,j,8);
        S310 = S(i,j,9);
        S120 = S(i,j,11);
        S220 = S(i,j,12);
        S030 = S(i,j,13);
        S130 = S(i,j,14);
        S040 = S(i,j,15);
        S101 = S(i,j,17);
        S201 = S(i,j,18);
        S301 = S(i,j,19);
        S102 = S(i,j,21);
        S202 = S(i,j,22);
        S003 = S(i,j,23);
        S103 = S(i,j,24);
        S004 = S(i,j,25);
        S011 = S(i,j,26);
        S111 = S(i,j,27);
        S211 = S(i,j,28);
        S021 = S(i,j,29);
        S121 = S(i,j,30);
        S031 = S(i,j,31);
        S012 = S(i,j,32);
        S112 = S(i,j,33);
        S013 = S(i,j,34);
        S022 = S(i,j,35);
        %
        D2 = delta2star3D(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
                          S101,S201,S301,S102,S202,S003,S103,S004,...
                          S011,S111,S211,S021,S121,S031,S012,S112,S013,S022);
        %dD2s = det(D2);
        dD2s = min([det(D2(1,1)) det(D2(1:2,1:2)) det(D2(1:3,1:3)) det(D2(1:4,1:4)) det(D2(1:5,1:5)) det(D2)]);
        if dD2s >= 0
            dDel2(i,j) = dD2s;
        else
            dDel2(i,j) = 0;
        end
    end
end
Z = dDel2;
contourf(X,Y,Z',50)
axis square
title('|\Delta_2^*|')
colorbar