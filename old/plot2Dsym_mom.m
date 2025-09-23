% plot moments vs. x at y = NP/2
nc = 3; nl = 3;
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

k=7;
subplot(nl,nc,k-1)
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
subplot(nl,nc,k-1)
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
subplot(nl,nc,k-1)
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

k=12;
subplot(nl,nc,k-3)
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

