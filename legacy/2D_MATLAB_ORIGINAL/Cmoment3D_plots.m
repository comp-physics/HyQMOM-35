% contour plots of various quantities 
%figure (10)
%cmin = 0; cmax = 100;

subplot(4,4,1)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,17);
contourf(X,Y,Z',50)
axis square
title('S_{101}')
colorbar

subplot(4,4,2)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,26);
contourf(X,Y,Z',50)
axis square
title('S_{011}')
colorbar

subplot(4,4,3)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,18);
contourf(X,Y,Z',50)
axis square
title('S_{201}')
colorbar

subplot(4,4,4)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,29);
contourf(X,Y,Z',50)
axis square
title('S_{021}')
colorbar

subplot(4,4,6)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,32);
contourf(X,Y,Z',50)
axis square
title('S_{012}')
colorbar

subplot(4,4,5)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,21);
contourf(X,Y,Z',50)
axis square
title('S_{102}')
colorbar

subplot(4,4,8)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,35)-1;
contourf(X,Y,Z',50)
axis square
title('S_{022}-1')
colorbar

subplot(4,4,7)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,22)-1;
contourf(X,Y,Z',50)
axis square
title('S_{202}-1')
colorbar

subplot(4,4,11)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,19);
contourf(X,Y,Z',50)
axis square
title('S_{301}')
colorbar

subplot(4,4,12)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,31);
contourf(X,Y,Z',50)
axis square
title('S_{031}')
colorbar

subplot(4,4,9)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,24);
contourf(X,Y,Z',50)
axis square
title('S_{103}')
colorbar

subplot(4,4,10)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,34);
contourf(X,Y,Z',50)
axis square
title('S_{013}')
colorbar

subplot(4,4,13)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,27);
contourf(X,Y,Z',50)
axis square
title('S_{111}')
colorbar

subplot(4,4,14)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,28);
contourf(X,Y,Z',50)
axis square
title('S_{211}')
colorbar

subplot(4,4,15)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,30);
contourf(X,Y,Z',50)
axis square
title('S_{121}')
colorbar

subplot(4,4,16)
[X,Y] = meshgrid(xm,ym);
Z = real(S(:,:,33)-S(:,:,7));
contourf(X,Y,Z',50)
axis square
title('S_{112}-S_{110}')
colorbar