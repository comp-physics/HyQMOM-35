% 3rd leading principal minor
clear
clc
close all

syms S110 S101 S011 S300 S210 S201 S120 S111 S102 S030 S021 S012 S003 S310 S301 S211 S130 S121 S112 S103 S031 S013 real
syms S400 S220 S202 S040 S022 S004 positive

S101=0;
S011=0;
S201=0;
S111=0;
S102=0;
S021=0;
S012=0;
S301=0;
S211=0;
S121=0;
S103=0;
S031=0;
S013=0;
S202=1;
S022=1;

dDel1 = -(S011^2 - 2*S011*S101*S110 + S101^2 + S110^2 - 1)

[E1,E2,E3,E4,E5,E6] = delta2star3D_permutation(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
                      S101,S201,S301,S102,S202,S003,S103,S004,S011,S111,S211,S021,S121,S031,S012,S112,S013,S022);
%

E6_3 = simplify(E6(1:3,1:3))

dE6_3 = det(E6_3)

solve(dE6_3,S112)


return

S110= 0.9;
S101= 0;
S011= 0;
[S110,S101,S011] = realizability_S2(S110,S101,S011)
%
S300=0;
S030=0;
S003=0;
S400=3;
S040=3;
S004=3;
%
S210=S110*S300;
S201=S101*S300;
S120=S110*S030;
S102=S101*S003;
S021=S011*S030;
S012=S011*S003;
%
A110 = ((S101-S011*S110)*S210+(S011-S101*S110)*S120)/(1-S110^2);
A101 = ((S110-S011*S101)*S201+(S011-S110*S101)*S102)/(1-S101^2);
A011 = ((S110-S101*S011)*S021+(S101-S110*S011)*S012)/(1-S011^2);
S111 = sign(A110)*min([abs(A110) abs(A101) abs(A011)]);
%
%
S310=S110*S400;
S301=S101*S400;
S211=0;
S130=S110*S040;
S121=0;
S112=0;
S103=S101*S004;
S031=S011*S040;
S013=S011*S004;

s220max = 1 + sqrt((S400-1)*(S040-1));
s202max = 1 + sqrt((S400-1)*(S004-1));
s022max = 1 + sqrt((S040-1)*(S004-1));

s220 = 0.9*s220max;
s202 = 0.9*s202max;
s022 = 0.1*s022max;

Del1 = [1 0    0    0;...
        0 1    S110 S101;...
        0 S110 1    S011;...
        0 S101 S011 1];

B = [1    S300 S210 S201;...
     S110 S210 S120 S111;...
     S101 S201 S111 S102;...
     1    S120 S030 S021;...
     S011 S111 S021 S012;...
     1    S102 S012 S003];

D = B/Del1*B';

d11 = D(1,1);
d22 = D(2,2);
d33 = D(3,3);
d44 = D(4,4);
d55 = D(5,5);
d66 = D(6,6);
d12 = D(1,2)
d13 = D(1,3);
d14 = D(1,4);
d16 = D(1,6);
d24 = D(2,4);
d36 = D(3,6);
d45 = D(4,5);
d46 = D(4,6);
d56 = D(5,6);

return

s220min = max([0 S110^2 1-sqrt((S400-1)*(S040-1)) d22+(S301-d13)^2/(S400-d11) d22+(S031-d45)^2/(S040-d44)]);
s202min = max([0 S101^2 1-sqrt((S400-1)*(S004-1)) d33+(S310-d12)^2/(S400-d11) d33+(S013-d56)^2/(S004-d66)]);
s022min = max([0 S011^2 1-sqrt((S040-1)*(S004-1)) d55+(S103-d36)^2/(S004-d66) d55+(S130-d24)^2/(S040-d44)]);

alpha = 0.5;

S220=s220min*alpha + (1-alpha)*s220max;
S202=s202min*alpha + (1-alpha)*s202max;
S022=s022min*alpha + (1-alpha)*s022max;

X0 = S220-d22;
Y0 = S202-d33;
Z0 = S022-d55;

x0 = s220-d22;
y0 = s202-d33;
z0 = s022-d55;

Xmax = s220max-d22;
Ymax = s202max-d33;
Zmax = s022max-d55;
Xmin = s220min-d22;
Ymin = s202min-d33;
Zmin = s022min-d55;

dDel1 = -(S011^2 - 2*S011*S101*S110 + S101^2 + S110^2 - 1)

[E1,~,~,~,~,~] = delta2star3D_permutation(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
                      S101,S201,S301,S102,S202,S003,S103,S004,S011,S111,S211,S021,S121,S031,S012,S112,S013,S022);
%
dE1 = det(E1)

[E0a,E0b,E0c,E0d,E0e,E0f] = delta2star3D_permutation(S300,S400,S110,S210,S310,S120,s220,S030,S130,S040,...
                      S101,S201,S301,S102,s202,S003,S103,S004,S011,S111,S211,S021,S121,S031,S012,S112,S013,s022);
%
dE04 = ([det(E0a(1:4,1:4)) det(E0b(1:4,1:4)) det(E0c(1:4,1:4)) det(E0d(1:4,1:4)) det(E0e(1:4,1:4)) det(E0f(1:4,1:4))])
dE05 = ([det(E0a(1:5,1:5)) det(E0b(1:5,1:5)) det(E0e(1:5,1:5))])
dE06 = det(E0a)

syms X Y positive

e11=E1(1,1);
e12=E1(1,2);
e13=E1(1,3);
e14=E1(1,4);
e15=E1(1,5);
e22=E1(2,2);
e23=E1(2,3);
e24=E1(2,4);
e25=E1(2,5);
e26=E1(2,6);
e33=E1(3,3);
e34=E1(3,4);
e35=E1(3,5);
e36=E1(3,6);
e44=E1(4,4);
e45=E1(4,5);
e55=E1(5,5);
e56=E1(5,6);
e66=E1(6,6);

ex = d22-d14;
ey = d33-d16;
ez = d55-d46;

N = 100;
x1 = linspace(Xmin,Xmax,N);
y1 = linspace(Ymin,Ymax,N);
[x,y] = meshgrid(x1,y1);
z1 = zeros(N,N);
z2 = z1;
z3 = z1;
for i =1:N
    for j=1:N
        X=x(i,j);
        Y=y(i,j);
        R = rootsZ_3D(X,Y,e11,e12,e13,e15,e23,e24,e25,e26,e34,e35,e36,e44,e45,e56,e66,ex,ey,ez);
        R = real(R);
        z1(i,j) = R(1);
        z2(i,j) = R(2);
        z3(i,j) = R(3);
        dD3 = detD3_3D(X,Y,R(1),e11,e12,e13,e15,e23,e24,e25,e26,e34,e35,e36,e44,e45,e56,e66,ex,ey,ez);
        dD31(i,j) = dD3;
        dD3 = detD3_3D(X,Y,R(2),e11,e12,e13,e15,e23,e24,e25,e26,e34,e35,e36,e44,e45,e56,e66,ex,ey,ez);
        dD32(i,j) = dD3;
        dD3 = detD3_3D(X,Y,R(3),e11,e12,e13,e15,e23,e24,e25,e26,e34,e35,e36,e44,e45,e56,e66,ex,ey,ez);
        dD33(i,j) = dD3;
        if max([abs(dD31(i,j)) abs(dD32(i,j)) abs(dD33(i,j))]) > 1.d-10 
            z1(i,j) = -1;
            z2(i,j) = -1;
            z3(i,j) = -1;
        end
    end
end

figure(1)
surface(x,y,z1)
view(3)
xlabel('X')
ylabel('Y')
zlabel('Z')
hold on
surface(x,y,z2)
surface(x,y,z3)
axis square
colormap sky
% xlim([Xmin Xmax])
% ylim([Ymin Ymax])
% zlim([Zmin Zmax])
xlim([0 Xmax])
ylim([0 Ymax])
zlim([0 Zmax])
plot3([Xmax Xmin],[Ymax Ymin],[Zmax Zmin],'ro-','MarkerFaceColor','r','LineWidth',1)
plot3([x0 X0],[y0 Y0],[z0 Z0],'bo-','MarkerFaceColor','b','LineWidth',1)
hold off

return

figure(2)
surface(x,y,dD31)
view(3)
xlabel('X')
ylabel('Y')
zlabel('Z')
hold on
surface(x,y,dD32)
surface(x,y,dD33)
axis square