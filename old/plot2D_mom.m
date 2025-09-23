% plot moments vs. x at y = NP/2
k=1;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{00}');
set(gca, 'FontSize', 18);
grid on;

k=2;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{10}');
set(gca, 'FontSize', 18);
grid on;

k=3;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{20}');
set(gca, 'FontSize', 18);
grid on;

k=4;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{30}');
set(gca, 'FontSize', 18);
grid on;

k=5;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{40}');
set(gca, 'FontSize', 18);
grid on;

k=6;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{01}');
set(gca, 'FontSize', 18);
grid on;

k=7;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{11}');
set(gca, 'FontSize', 18);
grid on;

k=8;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{21}');
set(gca, 'FontSize', 18);
grid on;

k=9;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{31}');
set(gca, 'FontSize', 18);
grid on;

k=10;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{02}');
set(gca, 'FontSize', 18);
grid on;

k=11;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{12}');
set(gca, 'FontSize', 18);
grid on;

k=12;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{22}');
set(gca, 'FontSize', 18);
grid on;

k=13;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{03}');
set(gca, 'FontSize', 18);
grid on;

k=14;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{13}');
set(gca, 'FontSize', 18);
grid on;

k=15;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{04}');
set(gca, 'FontSize', 18);
grid on;

k=16;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M5(i,i,1);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{50}');
set(gca, 'FontSize', 18);
grid on;

k=17;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M5(i,i,2);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{41}');
set(gca, 'FontSize', 18);
grid on;

k=18;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M5(i,i,3);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{32}');
set(gca, 'FontSize', 18);
grid on;

k=19;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M5(i,i,4);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{23}');
set(gca, 'FontSize', 18);
grid on;

k=20;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M5(i,i,5);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{14}');
set(gca, 'FontSize', 18);
grid on;