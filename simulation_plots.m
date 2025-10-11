function simulation_plots(plot_type, varargin)
% Comprehensive plotting functions for the 3D HyQMOM simulation
% All plotting functionality consolidated in one file
% 
% Usage:
%   simulation_plots('initial', xm, ym, M, C, S, M5, C5, S5, Np, enable_plots)
%   simulation_plots('time_evolution', S, C, xm, ym, Np, enable_plots)
%   simulation_plots('final', xm, ym, M, C, S, M5, C5, S5, Np, eig_data, enable_plots)
%   simulation_plots('final_time', xm, ym, M, Np, Nmom, Nmom5, enable_plots, txt)
%
% Note: If M, C, S arrays are 4D (with z-dimension), the middle z-slice is extracted automatically

switch lower(plot_type)
    case 'initial'
        plot_initial_conditions(varargin{:});
    case 'time_evolution'
        plot_time_evolution(varargin{:});
    case 'final'
        plot_final_results(varargin{:});
    case 'final_time'
        plot_final_time_standalone(varargin{:});
    otherwise
        error('Unknown plot type: %s. Valid types are: initial, time_evolution, final, final_time', plot_type);
end

end

% ========================================================================
% HELPER FUNCTIONS
% ========================================================================

function M_2d = extract_middle_z_slice(M)
% Extract middle z-slice from 4D array to get 3D array for plotting
% Input: M can be 3D [nx, ny, nmom] or 4D [nx, ny, nz, nmom]
% Output: M_2d is always 3D [nx, ny, nmom]
    if ndims(M) == 4
        nz = size(M, 3);
        k_mid = ceil(nz / 2);  % Middle z-index
        M_2d = squeeze(M(:, :, k_mid, :));
    else
        M_2d = M;  % Already 3D
    end
end

function plot_initial_conditions(xm, ym, M, C, S, M5, C5, S5, Np, enable_plots)
% Plot initial conditions for the 3D HyQMOM simulation
% This function creates figures 1-4 showing the initial state
if nargin < 10 || ~enable_plots
    return
end

% Extract middle z-slice if needed
M = extract_middle_z_slice(M);
C = extract_middle_z_slice(C);
S = extract_middle_z_slice(S);
M5 = extract_middle_z_slice(M5);
C5 = extract_middle_z_slice(C5);
S5 = extract_middle_z_slice(S5);

nmin = 1;
nmax = Np;
cc = 'k';

% Figure 1: Contour plots
figure(1)
contour_plots_3D(xm, ym, M, C, S, Np);
colormap sky

% Figure 2: Moment plots
figure(2)
plot_3Dsym_moments(xm, M, nmin, nmax, cc, Np);

% Figure 3: Central moment plots
figure(3)
plot_3Dsym_central(xm, M, C, C5, nmin, nmax, cc, Np);

% Figure 4: Standardized moment plots
figure(4)
plot_3Dsym_standardized(xm, S, S5, nmin, nmax, cc, Np);

pause(4)

end

function plot_time_evolution(S, C, xm, ym, Np, enable_plots)
% Plot time evolution during simulation (called in main loop)
% This function creates figure 10 showing current state during time stepping
if nargin < 6 || ~enable_plots
    return
end

% Extract middle z-slice if needed
S = extract_middle_z_slice(S);
C = extract_middle_z_slice(C);

figure(10)
Cmoment_plots_3D(xm, ym, S, Np);
colormap sky
pause(2)
hold off

end

function plot_final_results(xm, ym, M, C, S, M5, C5, S5, Np, eig_data, enable_plots)
% Plot final results after simulation completion
% This function creates figures 2-12 showing the final state and analysis
if nargin < 11 || ~enable_plots
    return
end

% Extract middle z-slice if needed
M = extract_middle_z_slice(M);
C = extract_middle_z_slice(C);
S = extract_middle_z_slice(S);
M5 = extract_middle_z_slice(M5);
C5 = extract_middle_z_slice(C5);
S5 = extract_middle_z_slice(S5);

nmin = 1;
nmax = Np;
cc = 'r';

% Figure 2: Final moment plots
figure(2)
clf;  % Clear this specific figure
plot_3Dsym_moments(xm, M, nmin, nmax, cc, Np);

% Figure 3: Final central moment plots
figure(3)
clf;
plot_3Dsym_central(xm, M, C, C5, nmin, nmax, cc, Np);

% Figure 4: Final standardized moment plots
figure(4)
clf;
plot_3Dsym_standardized(xm, S, S5, nmin, nmax, cc, Np);

% Figure 5: Eigenvalue comparison (x-direction)
figure(5)
clf;
Y1 = 0*xm;
Y2 = 0*xm;
Y3 = 0*xm;
Y4 = 0*xm;
for i=1:Np
    Y1(i) = eig_data.v6xmin(i,i);
    Y2(i) = eig_data.v6xmin(i,i);
    Y3(i) = eig_data.v6xmax(i,i);
    Y4(i) = eig_data.v6xmax(i,i);
end
plot(xm,Y1,'k',xm,Y2,'r',xm,Y3,'k',xm,Y4,'r')

% Figure 6: Eigenvalue plots (x-direction)
figure(6)
clf;
LAMXa = zeros(Np,6);
LAMXb = zeros(Np,6);
for i = 1:Np
    for kk = 1:6
        LAMXa(i,kk) = real(eig_data.lam6x(i,i,kk));
        LAMXb(i,kk) = real(eig_data.lam6z(i,i,kk));
    end
end
plot(xm,LAMXa,'o',xm,LAMXb,'p')

% Figure 7: Eigenvalue comparison (y-direction)
figure(7)
clf;
Y1 = 0*xm;
Y2 = 0*xm;
Y3 = 0*xm;
Y4 = 0*xm;
for i=1:Np
    Y1(i) = eig_data.v6ymin(i,i);
    Y2(i) = eig_data.v6ymin(i,i);
    Y3(i) = eig_data.v6ymax(i,i);
    Y4(i) = eig_data.v6ymax(i,i);
end
plot(xm,Y1,'k',xm,Y2,'r',xm,Y3,'k',xm,Y4,'r')

% Figure 8: Eigenvalue plots (y-direction)
figure(8)
clf;
LAMYa = zeros(Np,6);
LAMYb = zeros(Np,6);
for j = 1:Np
    for kk = 1:6
        LAMYa(j,kk) = real(eig_data.lam6y(j,j,kk));
        LAMYb(j,kk) = real(eig_data.lam6w(j,j,kk));
    end
end
plot(ym,LAMYa,'o',ym,LAMYb,'p')

% Figure 9: Final contour plots
figure(9)
clf;
contour_plots_3D(xm, ym, M, C, S, Np);
colormap sky

% Figure 10: Final C-moment plots
figure(10)
clf;
Cmoment_plots_3D(xm, ym, S, Np);
colormap sky

% Figure 11: Final S-moment plots
figure(11)
clf;
Smoment_plots_3D(xm, ym, S, Np);
colormap sky

% Figure 12: Hyperbolicity plots
figure(12)
clf;
hyperbolic_plots_3D(xm, ym, M, Np);
colormap sky

end

% ========================================================================
% INTERNAL PLOTTING FUNCTIONS
% ========================================================================

function contour_plots_3D(xm, ym, M, C, S, Np)
% Contour plots of various quantities (replaces contour3D_plots.m)
subplot(3,4,1)
[X,Y] = meshgrid(xm,ym);
Z = M(:,:,1);
contourf(X,Y,Z',50,'k')
axis square
title('Density')
colorbar

subplot(3,4,2)
[X,Y] = meshgrid(xm,ym);
Z = M(:,:,2)./M(:,:,1);
contourf(X,Y,Z',50,'k')
axis square
title('U velocity')
colorbar

subplot(3,4,3)
[X,Y] = meshgrid(xm,ym);
Z = M(:,:,6)./M(:,:,1);
contourf(X,Y,Z',50)
axis square
title('V velocity')
colorbar

subplot(3,4,4)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,7);
contourf(X,Y,Z',50)
axis square
title('S_{110}')
colorbar

subplot(3,4,5)
[X,Y] = meshgrid(xm,ym);
Z = C(:,:,3);
contourf(X,Y,Z',50)
axis square
title('C_{200}')
colorbar

subplot(3,4,6)
[X,Y] = meshgrid(xm,ym);
Z = C(:,:,10);
contourf(X,Y,Z',50)
axis square
title('C_{020}')
colorbar

subplot(3,4,7)
[X,Y] = meshgrid(xm,ym);
Z = C(:,:,20);
contourf(X,Y,Z',50)
axis square
title('C_{002}')
colorbar

subplot(3,4,8)
[X,Y] = meshgrid(xm,ym);
Z = 1 + 2*S(:,:,7).*S(:,:,17).*S(:,:,26) - S(:,:,7).^2 - S(:,:,17).^2 - S(:,:,26).^2;
contourf(X,Y,Z',50)
axis square
title('|\Delta_1|')
colorbar

subplot(3,4,9)
[X,Y] = meshgrid(xm,ym);
Z = max(0,S(:,:,5) - S(:,:,4).^2 - 1);
contourf(X,Y,Z',50)
axis square
title('H_{200}')
colorbar

subplot(3,4,10)
[X,Y] = meshgrid(xm,ym);
Z = max(0,S(:,:,15) - S(:,:,13).^2 - 1);
contourf(X,Y,Z',50)
axis square
title('H_{020}')
colorbar

subplot(3,4,11)
[X,Y] = meshgrid(xm,ym);
Z = max(0,S(:,:,25) - S(:,:,23).^2 - 1);
contourf(X,Y,Z',50)
axis square
title('H_{002}')
colorbar

subplot(3,4,12)
[X,Y] = meshgrid(xm,ym);
dDel2 = zeros(Np,Np);
for i = 1:Np
    for j = 1:Np
        % Use moment_idx for cleaner, self-documenting code
        S300 = S(i,j,moment_idx('M300')); S400 = S(i,j,moment_idx('M400')); 
        S110 = S(i,j,moment_idx('M110')); S210 = S(i,j,moment_idx('M210')); S310 = S(i,j,moment_idx('M310'));
        S120 = S(i,j,moment_idx('M120')); S220 = S(i,j,moment_idx('M220')); 
        S030 = S(i,j,moment_idx('M030')); S130 = S(i,j,moment_idx('M130')); S040 = S(i,j,moment_idx('M040'));
        S101 = S(i,j,moment_idx('M101')); S201 = S(i,j,moment_idx('M201')); S301 = S(i,j,moment_idx('M301')); 
        S102 = S(i,j,moment_idx('M102')); S202 = S(i,j,moment_idx('M202'));
        S003 = S(i,j,moment_idx('M003')); S103 = S(i,j,moment_idx('M103')); S004 = S(i,j,moment_idx('M004')); 
        S011 = S(i,j,moment_idx('M011')); S111 = S(i,j,moment_idx('M111'));
        S211 = S(i,j,moment_idx('M211')); S021 = S(i,j,moment_idx('M021')); S121 = S(i,j,moment_idx('M121')); 
        S031 = S(i,j,moment_idx('M031')); S012 = S(i,j,moment_idx('M012'));
        S112 = S(i,j,moment_idx('M112')); S013 = S(i,j,moment_idx('M013')); S022 = S(i,j,moment_idx('M022'));
        
        D2 = delta2star3D(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
                          S101,S201,S301,S102,S202,S003,S103,S004,...
                          S011,S111,S211,S021,S121,S031,S012,S112,S013,S022);
        dD2s = min([det(D2(1,1)) det(D2(1:2,1:2)) det(D2(1:3,1:3)) det(D2(1:4,1:4)) det(D2(1:5,1:5)) det(D2)]);
        dDel2(i,j) = max(0, dD2s);
    end
end
Z = dDel2;
contourf(X,Y,Z',50)
axis square
title('|\Delta_2^*|')
colorbar

end

function Cmoment_plots_3D(xm, ym, S, Np)
% C-moment plots (replaces Cmoment3D_plots.m)
subplot(4,4,1)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,17);
contourf(X,Y,Z',50)
axis square
title('S_{101}')
colorbar

subplot(4,4,2)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,26);
contourf(X,Y,Z',50)
axis square
title('S_{011}')
colorbar

subplot(4,4,3)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,18);
contourf(X,Y,Z',50)
axis square
title('S_{201}')
colorbar

subplot(4,4,4)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,29);
contourf(X,Y,Z',50)
axis square
title('S_{021}')
colorbar

subplot(4,4,5)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,21);
contourf(X,Y,Z',50)
axis square
title('S_{102}')
colorbar

subplot(4,4,6)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,32);
contourf(X,Y,Z',50)
axis square
title('S_{012}')
colorbar

subplot(4,4,7)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,22)-1;
contourf(X,Y,Z',50)
axis square
title('S_{202}-1')
colorbar

subplot(4,4,8)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,35)-1;
contourf(X,Y,Z',50)
axis square
title('S_{022}-1')
colorbar

subplot(4,4,9)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,24);
contourf(X,Y,Z',50)
axis square
title('S_{103}')
colorbar

subplot(4,4,10)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,34);
contourf(X,Y,Z',50)
axis square
title('S_{013}')
colorbar

subplot(4,4,11)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,19);
contourf(X,Y,Z',50)
axis square
title('S_{301}')
colorbar

subplot(4,4,12)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,31);
contourf(X,Y,Z',50)
axis square
title('S_{031}')
colorbar

subplot(4,4,13)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,27);
contourf(X,Y,Z',50)
axis square
title('S_{111}')
colorbar

subplot(4,4,14)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,28);
contourf(X,Y,Z',50)
axis square
title('S_{211}')
colorbar

subplot(4,4,15)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,30);
contourf(X,Y,Z',50)
axis square
title('S_{121}')
colorbar

subplot(4,4,16)
[X,Y] = meshgrid(xm,ym);
Z = real(S(:,:,33)-S(:,:,7));
contourf(X,Y,Z',50)
axis square
title('S_{112}-S_{110}')
colorbar

end

function Smoment_plots_3D(xm, ym, S, Np)
% S-moment plots (replaces Smoment3D_plots.m)
cmin = 0; cmax = 100;

subplot(3,4,1)
[X,Y] = meshgrid(xm,ym);
Z = max(-cmax,min(cmax,S(:,:,4)));
contourf(X,Y,Z',50)
axis square
title('S_{300}')
colorbar

subplot(3,4,2)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,8);
contourf(X,Y,Z',50)
axis square
title('S_{210}')
colorbar

subplot(3,4,3)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,11);
contourf(X,Y,Z',50)
axis square
title('S_{120}')
colorbar

subplot(3,4,4)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,32);
contourf(X,Y,Z',50)
axis square
title('S_{012}')
colorbar

subplot(3,4,5)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,13);
contourf(X,Y,Z',50)
axis square
title('S_{030}')
colorbar

subplot(3,4,6)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,5);
contourf(X,Y,Z',50)
axis square
title('S_{400}')
colorbar

subplot(3,4,7)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,9);
contourf(X,Y,Z',50)
axis square
title('S_{310}')
colorbar

subplot(3,4,8)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,12);
contourf(X,Y,Z',50)
axis square
title('S_{220}')
colorbar

subplot(3,4,9)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,14);
contourf(X,Y,Z',50)
axis square
title('S_{130}')
colorbar

subplot(3,4,10)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,25);
contourf(X,Y,Z',50)
axis square
title('S_{004}')
colorbar

subplot(3,4,11)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,33);
contourf(X,Y,Z',50)
axis square
title('S_{112}')
colorbar

subplot(3,4,12)
[X,Y] = meshgrid(xm,ym);
Z = S(:,:,35);
contourf(X,Y,Z',50)
axis square
title('S_{022}')
colorbar

end

function hyperbolic_plots_3D(xm, ym, M, Np)
% Hyperbolicity plots (replaces hyperbolic3D_plots.m)
% Use grid_eigenvalues to compute all eigenvalues efficiently
Nmom = size(M,3);
eig_data = grid_eigenvalues(M, Np, Nmom);

% Extract eigenvalue arrays for plotting
lam6x = eig_data.lam6x;  % UV plane
lam6y = eig_data.lam6y;  % VU plane  
lam6z = eig_data.lam6z;  % UW plane

% Plot eigenvalues for UV moments
[X,Y] = meshgrid(xm,ym);
for k = 1:6
    subplot(3,3,k)
    Z = lam6x(:,:,k);
    contourf(X,Y,Z',50)
    axis square
    title(['\lambda_' num2str(k)])
    colorbar
end

% Plot eigenvalue differences
subplot(3,3,7)
Z = lam6x(:,:,2) - lam6x(:,:,1);
contourf(X,Y,Z',50)
axis square
title('\lambda_2-\lambda_1')
colorbar

subplot(3,3,8)
Z = lam6x(:,:,4) - lam6x(:,:,3);
contourf(X,Y,Z',50)
axis square
title('\lambda_4-\lambda_3')
colorbar

subplot(3,3,9)
Z = lam6x(:,:,6) - lam6x(:,:,5);
contourf(X,Y,Z',50)
axis square
title('\lambda_6-\lambda_5')
colorbar

end

function plot_3Dsym_moments(xm, M, nmin, nmax, cc, Np)
% Plot moments vs. x at y = NP/2 (replaces plot3Dsym_mom.m)
% Uses 'hold on' to overlay initial (black) and final (red) results for comparison
% Note: M should be 3D [nx, ny, nmom] - z-slice already extracted
nc = 4; nl = 3;

% M000
subplot(nl,nc,1)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = M(i,i,1);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('M_{000}'); set(gca, 'FontSize', 18); grid on;

% M100
subplot(nl,nc,2)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = M(i,i,2);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('M_{100}'); set(gca, 'FontSize', 18); grid on;

% M200
subplot(nl,nc,3)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = M(i,i,3);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('M_{200}'); set(gca, 'FontSize', 18); grid on;

% M002
subplot(nl,nc,4)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = M(i,i,20);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('M_{002}'); set(gca, 'FontSize', 18); grid on;

% M300
subplot(nl,nc,5)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = M(i,i,4);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('M_{300}'); set(gca, 'FontSize', 18); grid on;

% M400
subplot(nl,nc,6)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = M(i,i,5);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('M_{400}'); set(gca, 'FontSize', 18); grid on;

% M110
subplot(nl,nc,7)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = M(i,i,7);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('M_{110}'); set(gca, 'FontSize', 18); grid on;

% M022
subplot(nl,nc,8)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = M(i,i,35);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('M_{022}'); set(gca, 'FontSize', 18); grid on;

% M210
subplot(nl,nc,9)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = M(i,i,8);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('M_{210}'); set(gca, 'FontSize', 18); grid on;

% M310
subplot(nl,nc,10)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = M(i,i,9);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('M_{310}'); set(gca, 'FontSize', 18); grid on;

% M220
subplot(nl,nc,11)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = M(i,i,12);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('M_{220}'); set(gca, 'FontSize', 18); grid on;

% M004
subplot(nl,nc,12)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = M(i,i,25);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('M_{004}'); set(gca, 'FontSize', 18); grid on;

end

function plot_3Dsym_central(xm, M, C, C5, nmin, nmax, cc, Np)
% Plot central moments vs. x at y = NP/2 (replaces plot3Dsym_C.m)
nl = 3; nc = 4;

% M000
subplot(nl,nc,1)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = M(i,i,1);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('M_{000}'); set(gca, 'FontSize', 18); grid on;

% u1 (mean velocity)
subplot(nl,nc,2)
hold on
umean = M(:,:,2)./M(:,:,1);
Y = zeros(size(xm));
for i=1:Np
    Y(i) = umean(i,i);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('u_1'); set(gca, 'FontSize', 18); grid on;

% C200
subplot(nl,nc,3)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = C(i,i,3);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('C_{200}'); set(gca, 'FontSize', 18); grid on;

% C300
subplot(nl,nc,4)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = C(i,i,4);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('C_{300}'); set(gca, 'FontSize', 18); grid on;

% C400
subplot(nl,nc,5)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = C(i,i,5);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('C_{400}'); set(gca, 'FontSize', 18); grid on;

% C110
subplot(nl,nc,6)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = C(i,i,7);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('C_{110}'); set(gca, 'FontSize', 18); grid on;

% C210
subplot(nl,nc,7)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = C(i,i,8);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('C_{210}'); set(gca, 'FontSize', 18); grid on;

% C310
subplot(nl,nc,8)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = C(i,i,9);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('C_{310}'); set(gca, 'FontSize', 18); grid on;

% C220
subplot(nl,nc,9)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = C(i,i,12);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('C_{220}'); set(gca, 'FontSize', 18); grid on;

% C500
subplot(nl,nc,10)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = C5(i,i,1);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('C_{500}'); set(gca, 'FontSize', 18); grid on;

% C410
subplot(nl,nc,11)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = C5(i,i,2);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('C_{410}'); set(gca, 'FontSize', 18); grid on;

% C320
subplot(nl,nc,12)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = C5(i,i,3);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('C_{320}'); set(gca, 'FontSize', 18); grid on;

end

function plot_3Dsym_standardized(xm, S, S5, nmin, nmax, cc, Np)
% Plot standardized moments vs. x at y = NP/2 (replaces plot3Dsym_S.m)
nl = 3; nc = 4;

% S300
subplot(nl,nc,1)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = S(i,i,4);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('S_{300}'); set(gca, 'FontSize', 18); grid on;

% S400
subplot(nl,nc,2)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = S(i,i,5);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('S_{400}'); set(gca, 'FontSize', 18); grid on;

% H400
subplot(nl,nc,3)
hold on
H40 = S(:,:,5) - S(:,:,4).^2 - 1;
Y = zeros(size(xm));
for i=1:Np
    Y(i) = H40(i,i);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('H_{400}'); set(gca, 'FontSize', 18); grid on;

% S500
subplot(nl,nc,4)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = S5(i,i,1);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('S_{500}'); set(gca, 'FontSize', 18); grid on;

% S110
subplot(nl,nc,5)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = S(i,i,7);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('S_{110}'); set(gca, 'FontSize', 18); grid on;

% S210
subplot(nl,nc,6)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = S(i,i,8);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('S_{210}'); set(gca, 'FontSize', 18); grid on;

% S310
subplot(nl,nc,7)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = S(i,i,9);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('S_{310}'); set(gca, 'FontSize', 18); grid on;

% S220
subplot(nl,nc,8)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = S(i,i,12);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('S_{220}'); set(gca, 'FontSize', 18); grid on;

% S410
subplot(nl,nc,9)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = S5(i,i,2);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('S_{410}'); set(gca, 'FontSize', 18); grid on;

% S320
subplot(nl,nc,10)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = S5(i,i,3);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('S_{320}'); set(gca, 'FontSize', 18); grid on;

% Delta1
subplot(nl,nc,11)
hold on
Y = zeros(size(xm));
for i=1:Np
    Y(i) = 1 + 2*S(i,i,7)*S(i,i,17)*S(i,i,26) - S(i,i,7)^2 - S(i,i,17)^2 - S(i,i,26)^2;
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('|\Delta_{1}|'); set(gca, 'FontSize', 18); grid on;

% Delta2*
subplot(nl,nc,12)
hold on
Y = zeros(size(xm));
for i = 1:Np
    S300 = S(i,i,4); S400 = S(i,i,5); S110 = S(i,i,7); S210 = S(i,i,8); S310 = S(i,i,9);
    S120 = S(i,i,11); S220 = S(i,i,12); S030 = S(i,i,13); S130 = S(i,i,14); S040 = S(i,i,15);
    S101 = S(i,i,17); S201 = S(i,i,18); S301 = S(i,i,19); S102 = S(i,i,21); S202 = S(i,i,22);
    S003 = S(i,i,23); S103 = S(i,i,24); S004 = S(i,i,25); S011 = S(i,i,26); S111 = S(i,i,27);
    S211 = S(i,i,28); S021 = S(i,i,29); S121 = S(i,i,30); S031 = S(i,i,31); S012 = S(i,i,32);
    S112 = S(i,i,33); S013 = S(i,i,34); S022 = S(i,i,35);
    
    D2 = delta2star3D(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
                      S101,S201,S301,S102,S202,S003,S103,S004,...
                      S011,S111,S211,S021,S121,S031,S012,S112,S013,S022);
    Y(i) = det(D2);
end
plot(xm(nmin:nmax),Y(nmin:nmax),'Color',cc,'LineWidth',2)
xlabel('x'); ylabel('|\Delta_{2}^*|'); set(gca, 'FontSize', 18); grid on;

end

function plot_final_time_standalone(xm, ym, M, Np, Nmom, Nmom5, enable_plots, txt)
% Standalone final time plotting (replaces plot_final_time.m)
% This function does postprocessing and creates final plots independently
if nargin < 7 || ~enable_plots
    return
end

% Clear and close existing plots
clc
close all

% Extract middle z-slice if needed
M = extract_middle_z_slice(M);

%% Postprocessing for plots
M5 = zeros(Np,Np,Nmom5);
C5 = zeros(Np,Np,Nmom5);
S5 = zeros(Np,Np,Nmom5);

for i = 1:Np
    for j = 1:Np
        MOM = zeros(Nmom,1);
        for k=1:Nmom
            MOM(k,1) = M(i,j,k);
        end
        [MM5,CC5,SS5] = Moments5_3D(MOM);
        for k=1:Nmom5
            M5(i,j,k) = MM5(k);
            C5(i,j,k) = CC5(k);
            S5(i,j,k) = SS5(k);
        end
    end
end

% Compute eigenvalues using grid_eigenvalues utility
eig_data = grid_eigenvalues(M, Np, Nmom);

% Extract eigenvalue arrays
lam6xa = eig_data.lam6x;   % UV plane (X-direction primary)
lam6xb = eig_data.lam6z;   % UW plane (X-direction secondary)
lam6ya = eig_data.lam6y;   % VU plane (Y-direction primary)
lam6yb = eig_data.lam6w;   % VW plane (Y-direction secondary)

% Save data if filename provided
if nargin >= 8 && ~isempty(txt)
    save(txt)
end

%% Create plots
nmin = 1;
nmax = Np;
cc = 'r';

% Compute C and S matrices for plotting
C = zeros(Np,Np,Nmom);
S = zeros(Np,Np,Nmom);
for i = 1:Np
    for j = 1:Np
        MOM = zeros(Nmom,1);
        for k=1:Nmom
            MOM(k,1) = M(i,j,k);
        end
        [CC,SS] = M2CS4_35(MOM);
        for k=1:Nmom
            C(i,j,k) = CC(k);
            S(i,j,k) = SS(k);
        end
    end
end

% Create velocity bounds (dummy values for compatibility)
v5xmin = zeros(Np,Np); v5xmax = zeros(Np,Np);
v5ymin = zeros(Np,Np); v5ymax = zeros(Np,Np);

% Figure 2: Moment plots
figure(2)
plot_3Dsym_moments(xm, M, nmin, nmax, cc, Np);

% Figure 3: Central moment plots
figure(3)
plot_3Dsym_central(xm, M, C, C5, nmin, nmax, cc, Np);

% Figure 4: Standardized moment plots
figure(4)
plot_3Dsym_standardized(xm, S, S5, nmin, nmax, cc, Np);

% Figure 5: Eigenvalue comparison (x-direction)
figure(5)
Y1 = 0*xm; Y2 = 0*xm; Y3 = 0*xm; Y4 = 0*xm;
for i=1:Np
    Y1(i) = v5xmin(i,i);
    Y2(i) = real(min(lam6xa(i,i,:)));
    Y3(i) = v5xmax(i,i);
    Y4(i) = real(max(lam6xa(i,i,:)));
end
plot(xm,Y1,'k',xm,Y2,'r',xm,Y3,'k',xm,Y4,'r')

% Figure 6: Eigenvalue plots (x-direction)
figure(6)
LAMXa = zeros(Np,6);
LAMXb = zeros(Np,6);
for i = 1:Np
    for k = 1:6
        LAMXa(i,k) = real(lam6xa(i,i,k));
        LAMXb(i,k) = real(lam6xb(i,i,k));
    end
end
plot(xm,LAMXa,'o',xm,LAMXb,'p')

% Figure 7: Eigenvalue comparison (y-direction)
figure(7)
Y1 = 0*xm; Y2 = 0*xm; Y3 = 0*xm; Y4 = 0*xm;
for i=1:Np
    Y1(i) = v5ymin(i,i);
    Y2(i) = real(min(lam6ya(i,i,:)));
    Y3(i) = v5ymax(i,i);
    Y4(i) = real(max(lam6ya(i,i,:)));
end
plot(xm,Y1,'k',xm,Y2,'r',xm,Y3,'k',xm,Y4,'r')

% Figure 8: Eigenvalue plots (y-direction)
figure(8)
LAMYa = zeros(Np,6);
LAMYb = zeros(Np,6);
for j = 1:Np
    for k = 1:6
        LAMYa(j,k) = real(lam6ya(j,j,k));
        LAMYb(j,k) = real(lam6yb(j,j,k));
    end
end
plot(ym,LAMYa,'o',ym,LAMYb,'p')

% Figure 9: Contour plots
figure(9)
contour_plots_3D(xm, ym, M, C, S, Np);
colormap sky

% Figure 10: C-moment plots
figure(10)
Cmoment_plots_3D(xm, ym, S, Np);
colormap sky

% Figure 11: S-moment plots
figure(11)
Smoment_plots_3D(xm, ym, S, Np);
colormap sky

% Figure 12: Hyperbolicity plots
figure(12)
hyperbolic_plots_3D(xm, ym, M, Np);
colormap sky

end