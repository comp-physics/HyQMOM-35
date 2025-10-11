% plot moments vs. x at y = NP/2
%figure(2)
nc = 4; nl = 3;
k=1;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,1);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{000}');
set(gca, 'FontSize', 18);
grid on;

k=2;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,2);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{100}');
set(gca, 'FontSize', 18);
grid on;

k=3;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,3);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{200}');
set(gca, 'FontSize', 18);
grid on;

k=4;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,20);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{002}');
set(gca, 'FontSize', 18);
grid on;

k=5;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,4);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{300}');
set(gca, 'FontSize', 18);
grid on;

k=6;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,5);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{400}');
set(gca, 'FontSize', 18);
grid on;

k=7;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,7);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{110}');
set(gca, 'FontSize', 18);
grid on;

k=8;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,35);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{022}');
set(gca, 'FontSize', 18);
grid on;

k=9;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,8);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{210}');
set(gca, 'FontSize', 18);
grid on;

k=10;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,9);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{310}');
set(gca, 'FontSize', 18);
grid on;

k=11;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,12);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{220}');
set(gca, 'FontSize', 18);
grid on;

k=12;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = M(i,i,25);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{004}');
set(gca, 'FontSize', 18);
grid on;