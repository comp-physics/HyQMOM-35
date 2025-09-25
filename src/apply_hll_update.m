function Mnp = apply_hll_update(M, F, vmin, vmax, dt, ds, direction)
% apply_hll_update Apply HLL flux update along specified direction
%
% Inputs:
%   M - moment array (Np x Np x Nmom)
%   F - flux array (Np x Np x Nmom)
%   vmin, vmax - velocity bounds (Np x Np)
%   dt - time step
%   ds - spatial step (dx or dy)
%   direction - 'x' or 'y'
%
% Output:
%   Mnp - updated moments (Np x Np x Nmom)

[Np, ~, Nmom] = size(M);
Mnp = M;

if strcmp(direction, 'x')
    % Update along x-direction (loop over j, process i-slices)
    parfor j = 1:Np
        VMIN = zeros(Np,1);
        VMAX = zeros(Np,1);
        MOM = zeros(Np,Nmom);
        FLUX = zeros(Np,Nmom);
        for i = 1:Np
            MOM(i,:) = squeeze(M(i,j,:));
            FLUX(i,:) = squeeze(F(i,j,:));
            VMIN(i,1) = vmin(i,j);
            VMAX(i,1) = vmax(i,j);
        end
        MNP = pas_HLL(MOM,FLUX,dt,ds,VMIN,VMAX);
        for i = 1:Np
            Mnp(i,j,:) = MNP(i,:);
        end
    end
elseif strcmp(direction, 'y')
    % Update along y-direction (loop over i, process j-slices)
    parfor i = 1:Np
        VMIN = zeros(Np,1);
        VMAX = zeros(Np,1);
        MOM = zeros(Np,Nmom);
        FLUX = zeros(Np,Nmom);
        for j = 1:Np
            MOM(j,:) = squeeze(M(i,j,:));
            FLUX(j,:) = squeeze(F(i,j,:));
            VMIN(j,1) = vmin(i,j);
            VMAX(j,1) = vmax(i,j);
        end
        MNP = pas_HLL(MOM,FLUX,dt,ds,VMIN,VMAX);
        for j = 1:Np
            Mnp(i,j,:) = MNP(j,:);
        end
    end
else
    error('Direction must be ''x'' or ''y''');
end

end
