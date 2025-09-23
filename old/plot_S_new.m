% plot standardized moments
k=1;
subplot(nl,nc,k)
hold on
plot(xm(nmin:nmax),M(nmin:nmax,k)','Color',cc,'LineWidth',2)
xlabel('x')
ylabel('M_{00}');
set(gca, 'FontSize', 18);
grid on;

k=2;
umean = M(:,2)./M(:,1);
subplot(nl,nc,k)
hold on
plot(xm(nmin:nmax),umean(nmin:nmax)','Color',cc,'LineWidth',2)
xlabel('x')
ylabel('u_1');
%ylim([-1 1])
set(gca, 'FontSize', 18);
grid on;

k=3;
H40 = S(:,5) - S(:,4).^2 - 1;
subplot(nl,nc,k)
hold on
plot(xm(nmin:nmax),H40(nmin:nmax)','Color',cc,'LineWidth',2)
xlabel('x')
ylabel('H_{40}');
%ylim([0 2])
set(gca, 'FontSize', 18);
grid on;

k=4;
subplot(nl,nc,k)
hold on
plot(xm(nmin:nmax),S(nmin:nmax,k)','Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{30}');
set(gca, 'FontSize', 18);
grid on;

k=5;
subplot(nl,nc,k)
hold on
plot(xm(nmin:nmax),S(nmin:nmax,k)','Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{40}');
set(gca, 'FontSize', 18);
grid on;

k=6;
vmean = M(:,6)./M(:,1);
subplot(nl,nc,k)
hold on
plot(xm(nmin:nmax),vmean(nmin:nmax)','Color',cc,'LineWidth',2)
xlabel('x')
ylabel('u_2');
%ylim([-1 1])
set(gca, 'FontSize', 18);
grid on;

k=7;
subplot(nl,nc,k)
hold on
plot(xm(nmin:nmax),S(nmin:nmax,k)','Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{11}');
%ylim([-1 1])
set(gca, 'FontSize', 18);
grid on;

k=8;
subplot(nl,nc,k)
hold on
plot(xm(nmin:nmax),S(nmin:nmax,k)','Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{21}');
set(gca, 'FontSize', 18);
grid on;

k=9;
subplot(nl,nc,k)
hold on
plot(xm(nmin:nmax),S(nmin:nmax,k)','Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{31}');
set(gca, 'FontSize', 18);
grid on;

k=10;
H40 = S(:,15) - S(:,13).^2 - 1;
subplot(nl,nc,k)
hold on
plot(xm(nmin:nmax),H40(nmin:nmax)','Color',cc,'LineWidth',2)
xlabel('x')
ylabel('H_{04}');
%ylim([0 2])
set(gca, 'FontSize', 18);
grid on;

k=11;
subplot(nl,nc,k)
hold on
plot(xm(nmin:nmax),S(nmin:nmax,k)','Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{12}');
%ylim([0 2])
set(gca, 'FontSize', 18);
grid on;

k=12;
subplot(nl,nc,k)
hold on
plot(xm(nmin:nmax),S(nmin:nmax,k)','Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{22}');
set(gca, 'FontSize', 18);
grid on;

k=13;
subplot(nl,nc,k)
hold on
plot(xm(nmin:nmax),S(nmin:nmax,k)','Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{03}');
set(gca, 'FontSize', 18);
grid on;

k=14;
subplot(nl,nc,k)
hold on
plot(xm(nmin:nmax),S(nmin:nmax,k)','Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{13}');
set(gca, 'FontSize', 18);
grid on;

k=15;
subplot(nl,nc,k)
hold on
plot(xm(nmin:nmax),S(nmin:nmax,k)','Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{04}');
set(gca, 'FontSize', 18);
grid on;

k=16;
subplot(nl,nc,k)
hold on
plot(xm(nmin:nmax),S5(nmin:nmax,1)','Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{50}');
set(gca, 'FontSize', 18);
grid on;

k=17;
subplot(nl,nc,k)
hold on
plot(xm(nmin:nmax),S5(nmin:nmax,2)','Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{41}');
set(gca, 'FontSize', 18);
grid on;

k=18;
subplot(nl,nc,k)
hold on
plot(xm(nmin:nmax),S5(nmin:nmax,3)','Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{32}');
set(gca, 'FontSize', 18);
grid on;

k=19;
subplot(nl,nc,k)
hold on
plot(xm(nmin:nmax),S5(nmin:nmax,4)','Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{23}');
set(gca, 'FontSize', 18);
grid on;

k=20;
subplot(nl,nc,k)
hold on
plot(xm(nmin:nmax),S5(nmin:nmax,5)','Color',cc,'LineWidth',2)
xlabel('x')
ylabel('S_{14}');
set(gca, 'FontSize', 18);
grid on;