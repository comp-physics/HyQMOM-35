% contour plots of various quantities 

cmin = 0; cmax = 100;

subplot(3,3,1)
[X,Y] = meshgrid(xm,ym);
Z = max(-cmax,min(cmax,S(:,:,4)));
contourf(X,Y,Z',50)
axis square
title('s_{30}')

subplot(3,3,2)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,8);
contourf(X,Y,Z',50)
axis square
title('s_{21}')

subplot(3,3,3)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,11);
contourf(X,Y,Z',50)
axis square
title('s_{12}')

subplot(3,3,4)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,13);
contourf(X,Y,Z',50)
axis square
title('s_{03}')

subplot(3,3,5)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,5);
contourf(X,Y,Z',50)
axis square
title('s_{40}')

subplot(3,3,6)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,9);
contourf(X,Y,Z',50)
axis square
title('s_{31}')

subplot(3,3,7)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,12);
contourf(X,Y,Z',50)
axis square
title('s_{22}')

subplot(3,3,8)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,14);
contourf(X,Y,Z',50)
axis square
title('s_{13}')

subplot(3,3,9)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,15);
contourf(X,Y,Z',50)
axis square
title('s_{04}')