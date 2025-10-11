% plot final figures
clc
close all

%% postprocessing for plots
for i = 1:Np
    for j = 1:Np
        MOM = zeros(Nmom,1);
        for k=1:Nmom
            MOM(k,1) = M(i,j,k);
        end
        [MM5,CC5,SS5] = Moments5_3D(MOM);
        for k=1:Nmom5
            M5(i,j,k) = MM5(k);
            C5(i,j,k) = CC5(k);
            S5(i,j,k) = SS5(k);
        end
    end
end

lam6xa = zeros(Np,Np,6);
lam6xb = zeros(Np,Np,6);
lam6ya = zeros(Np,Np,6);
lam6yb = zeros(Np,Np,6);
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
        J6 = jacobian6(m000,m010,m020,m030,m040,m100,m110,m120,m130,m200,m210,m220,m300,m310,m400);
        lam6xa(i,j,:) = eig(J6);
        J6 = jacobian6(m000,m001,m002,m003,m004,m100,m101,m102,m103,m200,m201,m202,m300,m301,m400);
        lam6xb(i,j,:) = eig(J6);
        J6 = jacobian6(m000,m100,m200,m300,m400,m010,m110,m210,m310,m020,m120,m220,m030,m130,m040);
        lam6ya(i,j,:) = eig(J6);
        J6 = jacobian6(m000,m001,m002,m003,m004,m010,m011,m012,m013,m020,m021,m022,m030,m031,m040);
        lam6yb(i,j,:) = eig(J6);
    end
end

% save date
save(txt)
%%

%% plots
nmin = 1;
nmax = Np;

cc = 'r';

figure(2)
%nc = 5; nl = 4;
plot3Dsym_mom

figure(3)
plot3Dsym_C

figure(4)
plot3Dsym_S

figure(5)
Y1 = 0*xm;
Y2 = 0*xm;
Y3 = 0*xm;
Y4 = 0*xm;
for i=1:Np
    Y1(i) = v5xmin(i,i);
    Y2(i) = v6xmin(i,i);
    Y3(i) = v5xmax(i,i);
    Y4(i) = v6xmax(i,i);
end
plot(xm,Y1,'k',xm,Y2,'r',xm,Y3,'k',xm,Y4,'r')

figure(6)
LAMXa = zeros(Np,6);
LAMXb = zeros(Np,6);
for i = 1:Np
    for k = 1:6
        LAMXa(i,k) = real(lam6xa(i,i,k));
        LAMXb(i,k) = real(lam6xb(i,i,k));
    end
end
plot(xm,LAMXa,'o',xm,LAMXb,'p')

figure(7)
Y1 = 0*xm;
Y2 = 0*xm;
Y3 = 0*xm;
Y4 = 0*xm;
for i=1:Np
    Y1(i) = v5ymin(i,i);
    Y2(i) = v6ymin(i,i);
    Y3(i) = v5ymax(i,i);
    Y4(i) = v6ymax(i,i);
end
plot(xm,Y1,'k',xm,Y2,'r',xm,Y3,'k',xm,Y4,'r')

figure(8)
LAMYa = zeros(Np,6);
LAMYb = zeros(Np,6);
for j = 1:Np
    for k = 1:6
        LAMYa(j,k) = real(lam6ya(j,j,k));
        LAMYb(j,k) = real(lam6yb(j,j,k));
    end
end
plot(ym,LAMYa,'o',ym,LAMYb,'p')

figure(9)
contour3D_plots
colormap sky

figure(10)
Cmoment3D_plots
colormap sky

figure(11)
Smoment3D_plots
colormap sky

figure(12)
hyperbolic3D_plots
colormap sky