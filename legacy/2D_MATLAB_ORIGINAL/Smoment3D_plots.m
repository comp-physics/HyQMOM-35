% contour plots of various quantities 
%figure (11)
cmin = 0; cmax = 100;

subplot(3,4,1)
[X,Y] = meshgrid(xm,ym);
Z = max(-cmax,min(cmax,S(:,:,4)));
contourf(X,Y,Z',50)
axis square
title('S_{300}')
colorbar

subplot(3,4,2)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,8);
contourf(X,Y,Z',50)
axis square
title('S_{210}')
colorbar

subplot(3,4,3)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,11);
contourf(X,Y,Z',50)
axis square
title('S_{120}')
colorbar

subplot(3,4,4)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,32);
contourf(X,Y,Z',50)
axis square
title('S_{012}')
colorbar

subplot(3,4,5)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,13);
contourf(X,Y,Z',50)
axis square
title('S_{030}')
colorbar

subplot(3,4,6)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,5);
contourf(X,Y,Z',50)
axis square
title('S_{400}')
colorbar

subplot(3,4,7)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,9);
contourf(X,Y,Z',50)
axis square
title('S_{310}')
colorbar

subplot(3,4,8)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,12);
contourf(X,Y,Z',50)
axis square
title('S_{220}')
colorbar

subplot(3,4,9)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,14);
contourf(X,Y,Z',50)
axis square
title('S_{130}')
colorbar

subplot(3,4,10)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,25);
contourf(X,Y,Z',50)
axis square
title('S_{004}')
colorbar

subplot(3,4,11)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,33);
contourf(X,Y,Z',50)
axis square
title('S_{112}')
colorbar

subplot(3,4,12)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,35);
contourf(X,Y,Z',50)
axis square
title('S_{022}')
colorbar