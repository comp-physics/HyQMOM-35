clear
clc

% syms umean vmean wmean real
%syms C200 C020 C002 positive
%syms C110 C101 C011 real
syms S300 S210 S201 S120 S111 S102 S030 S021 S012 S003 real
syms S400 S220 S202 S130 S040 S022 S004 positive
syms S310 S301 S211 S130 S121 S112 S103 S031 S013 real
% syms r110 r101 r011 real
% 
% assumeAlso(r110^2 < 1)
% assumeAlso(r101^2 < 1)
% assumeAlso(r011^2 < 1)
% assumeAlso(1 + 2*(r110*r101*r011) - r110^2 - r101^2 - r011^2 > 0)

syms C4 [5 5 5] real
syms A [3 3] real

% with rotation
%rhosr = [1 -1;1 1]*[sqrt(1+r) 0;0 sqrt(1-r)]*[1 1;-1 1]/2
% r110s = C110/sqrt(C200*C020);
% r101s = C101/sqrt(C200*C002);
% r011s = C011/sqrt(C020*C002);
% 
% rho = [1 r110 r101; r110 1 r011; r101 r011 1];
%rhosr = chol(rho,'real','nocheck');
%rhosr = sqrtm(rho);
%rhos = rhosr^2;

% r = r110;
% s = r101;
% t = r011;
% m = (r^2+s^2+t^2)/3;
% %if m > 0
% X = acos(-r*s*t/m/sqrt(m));
% lam1 = 1-2*sqrt(m)*cos(X);
% lam2 = 1-2*sqrt(m)*cos(X + 2*pi/3);
% lam3 = 1-2*sqrt(m)*cos(X + 4*pi/3);
% x1 = (s*(1-lam1)-r*t)/((1-lam1)^2-r^2);
% x2 = (s*(1-lam2)-r*t)/((1-lam2)^2-r^2);
% x3 = (s*(1-lam3)-r*t)/((1-lam3)^2-r^2);
% y1 = (t*(1-lam1)-r*s)/((1-lam1)^2-r^2);
% y2 = (t*(1-lam2)-r*s)/((1-lam2)^2-r^2);
% y3 = (t*(1-lam3)-r*s)/((1-lam3)^2-r^2);
% V1 =[x1, y1, -1]'/sqrt(x1^2+y1^2+1);
% V2 =[x2, y2, -1]'/sqrt(x2^2+y2^2+1);
% V3 =[x3, y3, -1]'/sqrt(x3^2+y3^2+1);
% Phi = [V1 V2 V3];
% rhosr = Phi*sqrt(diag([lam1,lam2,lam3]))*Phi';
% %
% 
% rhosr0 = double(subs(rhosr,[r110 r101 r011],[1 eps eps]))
% 
% return

%
C4(1,1,1) = 1;
%
C4(2,1,1) = 0;
C4(1,2,1) = 0;
C4(1,1,2) = 0;
%    
C4(3,1,1)=1;
C4(2,2,1)=0;
C4(2,1,2)=0;
C4(1,3,1)=1;
C4(1,2,2)=0;
C4(1,1,3)=1;
%
C4(4,1,1)=S300;
C4(3,2,1)=S210;
C4(3,1,2)=S201;
C4(2,3,1)=S120;
C4(2,2,2)=S111;
C4(2,1,3)=S102;
C4(1,4,1)=S030;
C4(1,3,2)=S021;
C4(1,2,3)=S012;
C4(1,1,4)=S003;
%
C4(5,1,1)=S400;
C4(4,2,1)=S310;
C4(4,1,2)=S301;
C4(3,3,1)=S220;
C4(3,2,2)=S211;
C4(3,1,3)=S202;
C4(2,4,1)=S130;
C4(2,3,2)=S121;
C4(2,2,3)=S112;
C4(2,1,4)=S103;
C4(1,5,1)=S040;
C4(1,4,2)=S031;
C4(1,3,3)=S022;
C4(1,2,4)=S013;
C4(1,1,5)=S004;

%simplify(rhosr*rhosr)
% B = [sqrt(C200) 0 0;0 sqrt(C020) 0; 0 0 sqrt(C002)]*rhosr;
% B = ...
% [    C200^(1/2),                                                  0,                                                                                   0; ...
% C020^(1/2)*r110,                      C020^(1/2)*(1 - r110^2)^(1/2),                                                                                   0; ...
% C002^(1/2)*r101, (C002^(1/2)*(r011 - r101*r110))/(1 - r110^2)^(1/2), C002^(1/2)*(((-r011^2 + 2*r011*r101*r110 - r101^2 - r110^2 + 1)/(1-r110^2))^(1/2)) ];
%B = subs(B,[r110 r101 r011],[r110s r101s r011s]);
%B = simplify(B)
% simplify(B*B')
 
%A = inv(B);

a11=A(1,1);
a12=A(1,2);
a13=A(1,3);
a21=A(2,1);
a22=A(2,2);
a23=A(2,3);
a31=A(3,1);
a32=A(3,2);
a33=A(3,3);
% b11=B(1,1);
% b12=B(1,2);
% b13=B(1,3);
% b21=B(2,1);
% b22=B(2,2);
% b23=B(2,3);
% b31=B(3,1);
% b32=B(3,2);
% b33=B(3,3);

S4 = 0*C4;

% find S4 from C4
for i=0:4 
    for j=0:4
        for k=0:4
            S4(i+1,j+1,k+1) = 0;
            if i+j+k <= 4
                for i1 = 0:i
                for i2 = 0:i-i1
                    i3 = i-i1-i2;
                    nchoosei1 = factorial(i)/factorial(i1)/factorial(i2)/factorial(i3);
                    for j1 = 0:j
                    for j2 = 0:j-j1
                        j3 = j-j1-j2;
                        nchoosej1 = factorial(j)/factorial(j1)/factorial(j2)/factorial(j3);
                        for k1 = 0:k
                        for k2 = 0:k-k1
                            k3 = k-k1-k2;
                            nchoosek1 = factorial(k)/factorial(k1)/factorial(k2)/factorial(k3);
                            S4(i+1,j+1,k+1) = S4(i+1,j+1,k+1) + nchoosei1*nchoosej1*nchoosek1*(a11^i1*a12^i2*a13^i3)*(a21^j1*a22^j2*a23^j3)*(a31^k1*a32^k2*a33^k3)*C4(1+i1+j1+k1,1+i2+j2+k2,1+i3+j3+k3);
                        end
                        end
                    end
                    end
                end
                end
            end
            simplify(S4(i+1,j+1,k+1));
        end   
    end
end

var = [A1_1,A1_2,A1_3,A2_1,A2_2,A2_3,A3_1,A3_2,A3_3,...
       S300,S210,S201,S120,S111,S102,S030,S021,S012,S003,...
       S400,S310,S301,S220,S211,S202,S130,S121,S112,S103,S040,S031,S022,S013,S004];

matlabFunction(S4,'File','S4toC4_3D_r','Vars',var)