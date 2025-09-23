%% fifth permutation E6: Z - Y
E = E6;

e11=E(1,1);
e12=E(1,2);
e13=E(1,3);
e14=E(1,4);
e15=E(1,5);
e16=E(1,6);
e22=E(2,2);
e23=E(2,3);
e24=E(2,4);
e25=E(2,5);
e26=E(2,6);
e33=E(3,3);
e34=E(3,4);
e35=E(3,5);
e36=E(3,6);
e44=E(4,4);
e45=E(4,5);
e46=E(4,6);
e55=E(5,5);
e56=E(5,6);
e66=E(6,6);

ex = -e22+e14;

x1 = linspace(Xmin,Xmax,N);
y1 = linspace(Ymin,Ymax,N);
C1 = ones(N,N,3);
C2 = C1;
C3 = C1;
[x,y] = meshgrid(x1,y1);
z1 = zeros(N,N);
z2 = z1;
z3 = z1;
z4 = z1;
z5 = z1;
z6 = z1;
for i =1:N
    for j=1:N
        Y = y(i,j);
        Z = x(i,j);
        R = rootsR_X_Y(Y,e11,e12,e13,e23,e24,e34,e44,ex);
        R = sort(real(R));
        z1(i,j) = R(1);
        z2(i,j) = R(2);
        z3(i,j) = R(3);
        if R(1) == R(2) || R(1) == R(3) || R(2) == R(3)
            C1(i,j,:) = [1 1 1];
            C2(i,j,:) = [1 1 1];
            C3(i,j,:) = [1 1 1];
        else
            C1(i,j,:) = [1 1 1];
            C2(i,j,:) = [1 0 0];
            C3(i,j,:) = [0 1 0];
        end
        R = rootsR_X_YZ(Y,Z,e11,e12,e13,e15,e23,e24,e25,e34,e35,e44,e45,ex);
        R = sort(real(R));
        z4(i,j) = R(1);
        z5(i,j) = R(2);
        z6(i,j) = R(3);
        if R(1) == R(2) || R(1) == R(3) || R(2) == R(3)
            C1(i,j,:) = [1 1 1];
            C2(i,j,:) = [1 1 1];
            C3(i,j,:) = [1 1 1];
        else
            C1(i,j,:) = [1 1 1];
            C2(i,j,:) = [1 0 0];
            C3(i,j,:) = [0 1 0];
        end
    end
end

%surface(x,y,z1,C1,'Linestyle','none')
surface(x,y,z2,C2,'Linestyle','none','FaceAlpha',FaceAlpha)
surface(x,y,z3,C3,'Linestyle','none','FaceAlpha',FaceAlpha)
%
% surface(x,y,z4,'FaceColor','y')
% surface(x,y,z5,'FaceColor','y')
% surface(x,y,z6,'FaceColor','y')