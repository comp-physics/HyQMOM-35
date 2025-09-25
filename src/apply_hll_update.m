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

% For y-direction, transpose to normalize to x-direction processing
flip = strcmp(direction, 'y');
if flip
    M = permute(M, [2 1 3]);
    F = permute(F, [2 1 3]);
    vmin = vmin.';
    vmax = vmax.';
end

Mnp = M;

% Process along first dimension (normalized to x-direction)
parfor j = 1:Np
    % Vectorized extraction
    MOM = squeeze(M(:,j,:));
    FLUX = squeeze(F(:,j,:));
    VMIN = vmin(:,j);
    VMAX = vmax(:,j);
    
    MNP = pas_HLL(MOM, FLUX, dt, ds, VMIN, VMAX);
    Mnp(:,j,:) = reshape(MNP, [size(MNP,1), 1, size(MNP,2)]);
end

% Transpose back if we processed y-direction
if flip
    Mnp = ipermute(Mnp, [2 1 3]);
end

end
