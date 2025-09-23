% plot standardized moments vs. x at y = NP/2
nl = 3; nc = 4;
k=4;
subplot(nl,nc,1)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = S(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{30}');
set(gca, 'FontSize', 18);
grid on;

k=5;
subplot(nl,nc,2)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = S(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{40}');
set(gca, 'FontSize', 18);
grid on;

k=3;
H40 = S(:,:,5) - S(:,:,4).^2 - 1;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = H40(i,i);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('H_{40}');
%ylim([0 2])
set(gca, 'FontSize', 18);
grid on;

k=4;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = S5(i,i,1);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{50}');
set(gca, 'FontSize', 18);
grid on;

k=7;
subplot(nl,nc,5)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = S(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{11}');
%ylim([-1 1])
set(gca, 'FontSize', 18);
grid on;

k=8;
subplot(nl,nc,6)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = S(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{21}');
set(gca, 'FontSize', 18);
grid on;

k=9;
subplot(nl,nc,7)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = S(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{31}');
set(gca, 'FontSize', 18);
grid on;

k=12;
subplot(nl,nc,8)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = S(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{22}');
set(gca, 'FontSize', 18);
grid on;

k=9;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = S5(i,i,2);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{41}');
set(gca, 'FontSize', 18);
grid on;

k=10;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = S5(i,i,3);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{32}');
set(gca, 'FontSize', 18);
grid on;

k=11;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = 1 - S(i,i,7).^2;
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('\Delta_{1}');
set(gca, 'FontSize', 18);
grid on;

k=12;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = S5(i,i,5);
end
for i = 1:Np
    s30 = S(i,i,4);
    s40 = S(i,i,5);
    s11 = S(i,i,7);
    s21 = S(i,i,8);
    s31 = S(i,i,9);
    s12 = S(i,i,11);
    s22 = S(i,i,12);
    s03 = S(i,i,13);
    s13 = S(i,i,14);
    s04 = S(i,i,15);
    Y(i) = delta2(s03,s04,s11,s12,s13,s21,s22,s30,s31,s40);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('\Delta_{2}');
set(gca, 'FontSize', 18);
grid on;