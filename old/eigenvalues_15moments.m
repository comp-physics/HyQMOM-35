% eigenvalues of 15 moment 2-D system from sec. 5.7

clear
clc
close all
% variables
syms m00 positive
syms m10 m01 m11 real
syms m20 m02 m22 positive
syms m30 m21 m12 m03 m31 m13 real
syms m40 m04 positive
% mean velocity
u1 = m10/m00;
u2 = m01/m00;
% 2nd-order central moments 
c20=m20/m00 - u1^2;
c02=m02/m00 - u2^2;
c11=m11/m00 - u1*u2;
%
assumeAlso(c20 > 0)
assumeAlso(c02 > 0)
assumeAlso(c11 > -sqrt(c20*c02))
assumeAlso(sqrt(c20*c02) > c11 )
%
assumeAlso(m40 - m30^2 - 1 > 0)
assumeAlso(m04 - m03^2 - 1 > 0)

% moments and closures for fluxes
Xx = [m00, m10, m20, m30, m40,   m01, m11, m21, m31,   m02, m12, m22,   m03, m13,   m04];
[Mx,My] = Flux_closure15(Xx);
% fluxes in the x direction
Fx = [m10, m20, m30, m40, Mx(1), m11, m21, m31, Mx(2), m12, m22, Mx(3), m13, Mx(4), Mx(5)];

% Jacobian matrices
Jx = jacobian(Fx,Xx);
Jx = simplify(Jx);
%spy(Jx)

matlabFunction(Jx(10:15,10:15),"File","jacobian6")

% change to standardized moments
syms s11 s30 s21 s12 s03 s31 s13 real
syms s40 s22 s04 h20 h02 positive
assumeAlso(s11 > -1)
assumeAlso(s11 < 1)
%
h20 = 0;
h02 = 0;
%s30 = 0; 
%s03 = 0;
%s40 = 1 + s30^2 + h20; 
%s04 = 1 + s03^2 + h02;
Jx0 = subs(Jx,[m00 m10 m01 m20 m11 m02 m30 m21 m12 m03 m40 m31 m22 m13 m04],[1 0 0 1 s11 1 s30 s21 s12 s03 s40 s31 1+s22 s13 s04]);
Jx0 = simplify(Jx0);
%spy(Jx0)

Jx1 = Jx0(1:5,1:5);
Jx2 = Jx0(6:9,6:9);
Jx3 = Jx0(10:15,10:15);

%Jx3a = subs(Jx3(1:3,1:3),s31,s11*s40+1.5*s30*(s21-s11*s30));
%Jx3b = subs(Jx3(4:6,4:6),[s31, s22, s13],[s11*s40+1.5*s30*(s21-s11*s30), 1, s11*s04+1.5*s03*(s12-s11*s03)]);

Jx3u = subs(Jx3,[s11 s12 s21 s31 s22 s13],[0 0 0 0 0 0]);
%spy(Jx3c)

Jx11 = diff(Jx3,s11);
%spy(Jx11)

Jx12 = diff(Jx3,s12);
%spy(Jx12)

Jx21 = diff(Jx3,s21);
%spy(Jx21)

Jx31 = diff(Jx3,s31);
%spy(Jx31)

Jx22 = diff(Jx3,s22);
%spy(Jx22)

Jx13 = diff(Jx3,s13);
%spy(Jx13)

%Jx3d = subs(Jx3b,[s11 s12 s21],[0 0 0]);
%spy(Jx3d)

% eigenvalues
%lambdax1  = eig(Jx1)
%lambdax2  = eig(Jx2)
%lambdax3  = eig(Jx3)

%return

syms x real

P5 = collect(det(x*eye(5)-Jx1),x)
P4 = collect(det(x*eye(4)-Jx2),x)
P6 = collect(det(x*eye(6)-Jx3),x)

%P3a = collect(det(x*eye(3)-Jx3a),x)
%P3b = collect(det(x*eye(3)-Jx3b),x)
P6u = collect(det(x*eye(6)-Jx3u),x)
%P3d = collect(det(x*eye(3)-Jx3d),x)


% define orthogonal polynomials
a00=0; a01=s03; a02=(a00+a01)/2;
b00=1; b01=1; b02=(s04-s03^2-1);

Q01 = x-a00;
Q02 = (x-a01)*Q01 - b01;
Q02 = collect(expand(Q02),x)

Q03 = (x-a02)*Q02 - b02*Q01;
Q03 = collect(expand(Q03),x)

R03 = (x-a02)*Q02 - 5/2*b02*Q01;
R03 = collect(expand(R03),x)

S03 = (x-a02)*Q02 - 3/2*b02*Q01;
S03 = collect(expand(S03),x)

% define orthogonal polynomials
a00=0; a10=s30; a20=(a00+a10)/2;
b00=1; b10=1; b20=(s40-s30^2-1);

Q10 = x-a00;
Q20 = (x-a10)*Q10 - b10;
Q20 = collect(expand(Q20),x)

Q30 = (x-a20)*Q20 - b20*Q10;
Q30 = collect(expand(Q30),x)

R30 = (x-a20)*Q20 - 5/2*b20*Q10;
R30 = collect(expand(R30),x)

Q40 = (x-a20)*Q30 - 3/2*b20*Q20;
Q40 = collect(expand(Q40),x)

[r,q] = polynomialReduce(P6,Q30,x);
r = simplify(r)
q = collect(simplify(q),x)

P6r = collect(simplify(P6-P6u),x)


return

% compare with P5, P3, P4
%dif5 =  collect(simplify(P5-Q20*R30),x)
%dif3 =  collect(simplify(P3-Q30),x)
%dif4 =  collect(simplify(P4-Q20*Q02),x)

% P412 = subs(P4,[s11 s12 s21 s40 s04],[1 s03 s30 1+s30^2 1+s03^2]);
% P412 = collect(simplify(P412),x) 
% 
% solve(P412,x)
syms h positive

P412 = subs(P4,[s30 s21 s12 s03],[0 0 0 0]);
P412 = collect(simplify(P412),x) 

%x0 = solve(P412,x)
S2 = x^2

[r,q] = polynomialReduce(P412,S2,x);
r = simplify(r)
q = collect(simplify(q),x)

% P412 = subs(P4,[s12 s04],[s11*s03 1+s03^2]);
% P412 = collect(simplify(P412),x) 
% 
% solve(P412,x)

return

% syms s50 real
% f = expand(P412);
% f = collect(f,x);
% f = subs(f,[x^5 x^4 x^3 x^2 x],[s50 s40 s30 1 0]);
% f1 = simplify(f)
% 
% f = expand(x*P412);
% f = collect(f,x);
% f = subs(f,[x^5 x^4 x^3 x^2 x],[s50 s40 s30 1 0]);
% f2 = simplify(f)
% s50s = solve(f2,s50)

P4b = subs(P4,[s12],[s11*s03]);
P4b = collect(simplify(P4b),x) 

S2 = x^2 + (3*s03*s11 - 4*s12)*x + (9*s03^2*s11^2)/4 - 6*s03*s11*s12 - s11^2 + (15*s12^2)/4

[r,q] = polynomialReduce(P4,S2,x)
r = simplify(r)
q = collect(simplify(q),x)

% syms b3 real

% P4 = x^4 + ((3*s03*s11)/2 - s30 - b3*s12 - (5*s12)/2 + b3*s03*s11)*x^3 + ((9*s03^2*s11^2)/4 - (s03*s21)/2 + (5*s12*s30)/2 - b3*s11^2 + (3*b3*s12^2)/2 - (3*s04*s11^2)/2 + s11^2/2 + (3*s12^2)/2 + b3*s04*s11^2 + b3*s12*s30 - (9*s03*s11*s12)/4 - s03*s11*s30 - (5*b3*s03*s11*s12)/2 - b3*s03*s11*s30 - 1)*x^2 + ((5*s12)/2 + b3*s12 - (3*s03*s11)/2 - (s11^2*s30)/2 - (3*s12^2*s30)/2 - (3*b3*s12^2*s30)/2 + (3*s04*s11^2*s30)/2 - (9*s03^2*s11^2*s30)/4 - b3*s03*s11 + b3*s11*s21 + (s03*s12*s21)/2 + (b3*s03*s12*s21)/2 - b3*s04*s11*s21 + (7*s03*s11*s12*s30)/4 + (b3*s03^2*s11*s21)/2 - (b3*s03^2*s11^2*s30)/2 + 2*b3*s03*s11*s12*s30)*x + b3*s11^2 - (9*s03^2*s11^2)/4 - (3*b3*s12^2)/2 + (3*s04*s11^2)/2 - s11^2/2 - (3*s12^2)/2 - b3*s04*s11^2 + (s03*s11^2*s21)/2 - (s03*s11^3*s30)/2 + (9*s03*s11*s12)/4 + (5*b3*s03*s11*s12)/2 + b3*s11*s12*s21 - b3*s03*s11^2*s21 + b3*s03*s11^3*s30 - b3*s11^2*s12*s30;
% P4p= diff(P4,x) 
% 
% %P4m = s04 + ((3*s03*s11)/2 - s30 - b3*s12 - (5*s12)/2 + b3*s03*s11)*s03 + ((9*s03^2*s11^2)/4 - (s03*s21)/2 + (5*s12*s30)/2 - b3*s11^2 + (3*b3*s12^2)/2 - (3*s04*s11^2)/2 + s11^2/2 + (3*s12^2)/2 + b3*s04*s11^2 + b3*s12*s30 - (9*s03*s11*s12)/4 - s03*s11*s30 - (5*b3*s03*s11*s12)/2 - b3*s03*s11*s30 - 1) + ((5*s12)/2 + b3*s12 - (3*s03*s11)/2 - (s11^2*s30)/2 - (3*s12^2*s30)/2 - (3*b3*s12^2*s30)/2 + (3*s04*s11^2*s30)/2 - (9*s03^2*s11^2*s30)/4 - b3*s03*s11 + b3*s11*s21 + (s03*s12*s21)/2 + (b3*s03*s12*s21)/2 - b3*s04*s11*s21 + (7*s03*s11*s12*s30)/4 + (b3*s03^2*s11*s21)/2 - (b3*s03^2*s11^2*s30)/2 + 2*b3*s03*s11*s12*s30)*0 + b3*s11^2 - (9*s03^2*s11^2)/4 - (3*b3*s12^2)/2 + (3*s04*s11^2)/2 - s11^2/2 - (3*s12^2)/2 - b3*s04*s11^2 + (s03*s11^2*s21)/2 - (s03*s11^3*s30)/2 + (9*s03*s11*s12)/4 + (5*b3*s03*s11*s12)/2 + b3*s11*s12*s21 - b3*s03*s11^2*s21 + b3*s03*s11^3*s30 - b3*s11^2*s12*s30;
% 
% P4p = collect(simplify(P4p),x)
% 
% %b3s = solve(P4m,b3)
% 
% %P4px = 4*s04 + ((3*s03*s11)/2 - s30 - b3*s12 - (5*s12)/2 + b3*s03*s11)*3*s03 + ((9*s03^2*s11^2)/4 - (s03*s21)/2 + (5*s12*s30)/2 - b3*s11^2 + (3*b3*s12^2)/2 - (3*s04*s11^2)/2 + s11^2/2 + (3*s12^2)/2 + b3*s04*s11^2 + b3*s12*s30 - (9*s03*s11*s12)/4 - s03*s11*s30 - (5*b3*s03*s11*s12)/2 - b3*s03*s11*s30 - 1)*2 + ((5*s12)/2 + b3*s12 - (3*s03*s11)/2 - (s11^2*s30)/2 - (3*s12^2*s30)/2 - (3*b3*s12^2*s30)/2 + (3*s04*s11^2*s30)/2 - (9*s03^2*s11^2*s30)/4 - b3*s03*s11 + b3*s11*s21 + (s03*s12*s21)/2 + (b3*s03*s12*s21)/2 - b3*s04*s11*s21 + (7*s03*s11*s12*s30)/4 + (b3*s03^2*s11*s21)/2 - (b3*s03^2*s11^2*s30)/2 + 2*b3*s03*s11*s12*s30);%*x + b3*s11^2 - (9*s03^2*s11^2)/4 - (3*b3*s12^2)/2 + (3*s04*s11^2)/2 - s11^2/2 - (3*s12^2)/2 - b3*s04*s11^2 + (s03*s11^2*s21)/2 - (s03*s11^3*s30)/2 + (9*s03*s11*s12)/4 + (5*b3*s03*s11*s12)/2 + b3*s11*s12*s21 - b3*s03*s11^2*s21 + b3*s03*s11^3*s30 - b3*s11^2*s12*s30;
% 
% 
% %P4p = (5*s12)/2 + b3*s12 - (3*s03*s11)/2 - 2*x*((s03*s21)/2 - (9*s03^2*s11^2)/4 - (5*s12*s30)/2 + b3*s11^2 - (3*b3*s12^2)/2 + (3*s04*s11^2)/2 - s11^2/2 - (3*s12^2)/2 - b3*s04*s11^2 - b3*s12*s30 + (9*s03*s11*s12)/4 + s03*s11*s30 + (5*b3*s03*s11*s12)/2 + b3*s03*s11*s30 + 1) - (s11^2*s30)/2 - (3*s12^2*s30)/2 + 4*x^3 - 3*x^2*((5*s12)/2 + s30 + b3*s12 - (3*s03*s11)/2 - b3*s03*s11) - (3*b3*s12^2*s30)/2 + (3*s04*s11^2*s30)/2 - (9*s03^2*s11^2*s30)/4 - b3*s03*s11 + b3*s11*s21 + (s03*s12*s21)/2 + (b3*s03*s12*s21)/2 - b3*s04*s11*s21 + (7*s03*s11*s12*s30)/4 + (b3*s03^2*s11*s21)/2 - (b3*s03^2*s11^2*s30)/2 + 2*b3*s03*s11*s12*s30
% 
% 
% P4p = 4*s03 + ((9*s03*s11)/2 - 3*s30 - 3*b3*s12 - (15*s12)/2 + 3*b3*s03*s11) + ((9*s03^2*s11^2)/2 - s03*s21 + 5*s12*s30 - 2*b3*s11^2 + 3*b3*s12^2 - 3*s04*s11^2 + s11^2 + 3*s12^2 + 2*b3*s04*s11^2 + 2*b3*s12*s30 - (9*s03*s11*s12)/2 - 2*s03*s11*s30 - 5*b3*s03*s11*s12 - 2*b3*s03*s11*s30 - 2)*0 + (5*s12)/2 + b3*s12 - (3*s03*s11)/2 - (s11^2*s30)/2 - (3*s12^2*s30)/2 - (3*b3*s12^2*s30)/2 + (3*s04*s11^2*s30)/2 - (9*s03^2*s11^2*s30)/4 - b3*s03*s11 + b3*s11*s21 + (s03*s12*s21)/2 + (b3*s03*s12*s21)/2 - b3*s04*s11*s21 + (7*s03*s11*s12*s30)/4 + (b3*s03^2*s11*s21)/2 - (b3*s03^2*s11^2*s30)/2 + 2*b3*s03*s11*s12*s30;
% P4px = 4*s04 + ((9*s03*s11)/2 - 3*s30 - 3*b3*s12 - (15*s12)/2 + 3*b3*s03*s11)*s03 + ((9*s03^2*s11^2)/2 - s03*s21 + 5*s12*s30 - 2*b3*s11^2 + 3*b3*s12^2 - 3*s04*s11^2 + s11^2 + 3*s12^2 + 2*b3*s04*s11^2 + 2*b3*s12*s30 - (9*s03*s11*s12)/2 - 2*s03*s11*s30 - 5*b3*s03*s11*s12 - 2*b3*s03*s11*s30 - 2);% + (5*s12)/2 + b3*s12 - (3*s03*s11)/2 - (s11^2*s30)/2 - (3*s12^2*s30)/2 - (3*b3*s12^2*s30)/2 + (3*s04*s11^2*s30)/2 - (9*s03^2*s11^2*s30)/4 - b3*s03*s11 + b3*s11*s21 + (s03*s12*s21)/2 + (b3*s03*s12*s21)/2 - b3*s04*s11*s21 + (7*s03*s11*s12*s30)/4 + (b3*s03^2*s11*s21)/2 - (b3*s03^2*s11^2*s30)/2 + 2*b3*s03*s11*s12*s30
% 
% P4p = simplify(P4p)
% P4px = simplify(P4px)
% 
% b3p = solve(P4p,b3)
% b3x = solve(P4px,b3)

% P4s = -((- x^2 + s30*x + 1)*(9*s03^2*s11^2 - 24*s03*s11*s12 + 12*s03*s11*x - 4*s11^2 + 15*s12^2 - 16*s12*x + 4*x^2))/4;
% S2s = 9*s03^2*s11^2 - 24*s03*s11*s12 + 12*s03*s11*x - 4*s11^2 + 15*s12^2 - 16*s12*x + 4*x^2;
% 
% solve(S2s,x)
 