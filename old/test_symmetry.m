function [Diff,MaxDiff] = test_symmetry(M,Np)
% Check symmetry of moment matrix
%
Diff = zeros(Np,5);
for i = 1:Np
    Diff(i,1) = M(i,i,1) - M(Np+1-i,Np+1-i,1);
    Diff(i,2) = M(i,i,2) + M(Np+1-i,Np+1-i,2);
    Diff(i,3) = M(i,i,3) - M(Np+1-i,Np+1-i,3);
    Diff(i,4) = M(i,i,4) + M(Np+1-i,Np+1-i,4);
    Diff(i,5) = M(i,i,5) - M(Np+1-i,Np+1-i,5);
end
MaxDiff = zeros(5,1);
for k = 1:5
    MaxDiff(k) = max(Diff(:,k));
end
end