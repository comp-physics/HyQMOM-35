function [lam6xa, lam6xb, lam6ya, lam6yb] = compute_jacobian_eigenvalues(M)
% compute_jacobian_eigenvalues Compute 6x6 Jacobian eigenvalues for all grid points
%
% Input:
%   M - moment array (Np x Np x Nmom)
%
% Outputs:
%   lam6xa, lam6xb, lam6ya, lam6yb - eigenvalue arrays (Np x Np x 6)

[Np, ~, ~] = size(M);
lam6xa = zeros(Np,Np,6);
lam6xb = zeros(Np,Np,6);
lam6ya = zeros(Np,Np,6);
lam6yb = zeros(Np,Np,6);

for i = 1:Np
    for j = 1:Np
        M1 = squeeze(M(i,j,:));
        
        % Use template-based Jacobian computation
        x_templates = jacobian_templates('x');
        y_templates = jacobian_templates('y');
        
        % X-direction eigenvalues
        args_xa = num2cell(M1(x_templates{1}));
        J6 = jacobian6(args_xa{:});
        lam6xa(i,j,:) = eig(J6);
        
        args_xb = num2cell(M1(x_templates{2}));
        J6 = jacobian6(args_xb{:});
        lam6xb(i,j,:) = eig(J6);
        
        % Y-direction eigenvalues
        args_ya = num2cell(M1(y_templates{1}));
        J6 = jacobian6(args_ya{:});
        lam6ya(i,j,:) = eig(J6);
        
        args_yb = num2cell(M1(y_templates{2}));
        J6 = jacobian6(args_yb{:});
        lam6yb(i,j,:) = eig(J6);
    end
end

end
