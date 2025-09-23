function [Diff,MaxDiff] = test_symmetry_1Dx(M,Np,j)
% Check symmetry of moment matrix in x direction
%
Diff = zeros(Np,5);
for i = 1:Np
    Diff(i,1) = M(i,j,1) - M(Np+1-i,j,1);
    Diff(i,2) = M(i,j,2) + M(Np+1-i,j,2);
    Diff(i,3) = M(i,j,3) - M(Np+1-i,j,3);
    Diff(i,4) = M(i,j,4) + M(Np+1-i,j,4);
    Diff(i,5) = M(i,j,5) - M(Np+1-i,j,5);
end
MaxDiff = zeros(5,1);
for k = 1:5
    MaxDiff(k) = max(Diff(:,k));
end
end