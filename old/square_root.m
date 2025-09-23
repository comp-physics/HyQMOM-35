clc
clear

syms e11 e44 X Y positive
syms e12 e13 e23 e24 e34 ex real

v1 = [X+ex e24 e34]';

C = inv([e11 e12 e13; e12 X e23; e13 e23 Y]);

D = sqrtm(C);

v2 = D*v1;

mv2 = e44 - v2'*v2;
mV2 = collect(mv2,X)

