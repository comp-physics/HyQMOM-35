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
ylabel('S_{300}');
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
ylabel('S_{400}');
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
ylabel('H_{400}');
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
ylabel('S_{500}');
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
ylabel('S_{110}');
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
ylabel('S_{210}');
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
ylabel('S_{310}');
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
ylabel('S_{220}');
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
ylabel('S_{410}');
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
ylabel('S_{320}');
set(gca, 'FontSize', 18);
grid on;

k=11;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i=1:Np
    Y(i) = 1 + 2*S(i,i,7)*S(i,i,17)*S(i,i,26) - S(i,i,7)^2 - S(i,i,17)^2 - S(i,i,26)^2;
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('|\Delta_{1}|');
set(gca, 'FontSize', 18);
grid on;

k=12;
subplot(nl,nc,k)
hold on
Y = 0*xm;
for i = 1:Np
    S300 = S(i,i,4);
    S400 = S(i,i,5);
    S110 = S(i,i,7);
    S210 = S(i,i,8);
    S310 = S(i,i,9);
    S120 = S(i,i,11);
    S220 = S(i,i,12);
    S030 = S(i,i,13);
    S130 = S(i,i,14);
    S040 = S(i,i,15);
    S101 = S(i,i,17);
    S201 = S(i,i,18);
    S301 = S(i,i,19);
    S102 = S(i,i,21);
    S202 = S(i,i,22);
    S003 = S(i,i,23);
    S103 = S(i,i,24);
    S004 = S(i,i,25);
    S011 = S(i,i,26);
    S111 = S(i,i,27);
    S211 = S(i,i,28);
    S021 = S(i,i,29);
    S121 = S(i,i,30);
    S031 = S(i,i,31);
    S012 = S(i,i,32);
    S112 = S(i,i,33);
    S013 = S(i,i,34);
    S022 = S(i,i,35);
    %
    D2 = delta2star3D(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
                      S101,S201,S301,S102,S202,S003,S103,S004,...
                      S011,S111,S211,S021,S121,S031,S012,S112,S013,S022);
    Y(i)  = det(D2);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x')
ylabel('|\Delta_{2}^*|');
set(gca, 'FontSize', 18);
grid on;