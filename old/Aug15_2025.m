clear
clc

% syms x y z real
% 
% x = y*z + sqrt((1-y^2)*(1-z^2));
% 
% A = [1 x y; x 1 z; y z 1];
% 
% [V,D] = eig(A)
% 
% X = (y*((y^2 - 1)*(z^2 - 1))^(1/2) - z + y^2*z)*V(:,1)
% 
% simplify(A*X)

syms a b c positive

% A=sym(zeros(6,6)); 
% A(1,1)=b; A(1,2)=c; 
% A(2,1)=a; A(2,3)=b; A(3,2)=a;
% A(3,4)=c; A(4,3)=a; A(4,5)=c; A(5,4)=a; A(5,6)=b; A(6,5)=b; A(6,6)=c;
% A
% Ap=A(2:6,2:6)
% rank(Ap)
% 
% [L,U,P] = lu(Ap)

B = [a^2 -c*a*b; -c*a*b b^2];
Bs = sqrtm(B)

[V,D]=eig(B)

syms x1 x2 real
X = [x1; x2];

LHS = simplify(X'*B*X)

