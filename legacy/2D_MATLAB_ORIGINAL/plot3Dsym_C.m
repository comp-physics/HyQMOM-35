% plot central moments vs. x at y = NP/2
nl = 3; nc = 4;
k=1;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{000}');
set(gca, 'FontSize', 18);
grid on;

k=2;
umean = M(:,:,2)./M(:,:,1);
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = umean(i,i);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('u_1');
set(gca, 'FontSize', 18);
grid on;

k=3;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = C(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('C_{200}');
set(gca, 'FontSize', 18);
grid on;

k=4;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = C(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('C_{300}');
set(gca, 'FontSize', 18);
grid on;

k=5;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = C(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('C_{400}');
set(gca, 'FontSize', 18);
grid on;

k=7;
subplot(nl,nc,k-1)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = C(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('C_{110}');
set(gca, 'FontSize', 18);
grid on;

k=8;
subplot(nl,nc,k-1)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = C(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('C_{210}');
set(gca, 'FontSize', 18);
grid on;

k=9;
subplot(nl,nc,k-1)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = C(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('C_{310}');
set(gca, 'FontSize', 18);
grid on;

k=12;
subplot(nl,nc,k-3)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = C(i,i,k);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('C_{220}');
set(gca, 'FontSize', 18);
grid on;

k=10;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = C5(i,i,1);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('C_{500}');
set(gca, 'FontSize', 18);
grid on;

k=11;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = C5(i,i,2);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('C_{410}');
set(gca, 'FontSize', 18);
grid on;

k=12;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = C5(i,i,3);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('C_{320}');
set(gca, 'FontSize', 18);
grid on;
