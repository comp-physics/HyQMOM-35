% 2-D Riemann solver for 3-D HyQMOM using HLL and explicit Euler
% code restricted to N=4 in 3-D with 35 moments
%
% M = [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,...
%      M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,...
%      M031,M012,M112,M013,M022]
%
% initial conditions can be changed to test different scenario
%
% This version limits the energy fluxes and checks realizability and
% hyperbolicity
%
clc
clear 
close all

% Knudsen number (>= 0.001 to avoid long simulations)
Kn = 1/1;

% Mach number (for impinging jets with velocity u and temperature Theta)
Ma = 0;  % (= u/sqrt(Theta))

% final time (<= 0.075 to keep waves in 2-D box - smaller for large Ma)
tmax = 0.1 ;

% flag for 2-D case (use only if S101=S011=0) if flag2D == 1
flag2D = 0;

%% 2-D space discretization: square domain
Np = 40;
CFL = 0.5;
xmin = -0.5;
xmax = 0.5;
ymin = -0.5;
ymax = 0.5;
x = xmin + (xmax-xmin)*linspace(0,1,Np+1)';
y = ymin + (ymax-ymin)*linspace(0,1,Np+1)';
dx = (xmax-xmin)/Np;
xm = x(1:Np)+dx/2;
dy = (ymax-ymin)/Np;
ym = y(1:Np)+dy/2;
%%

% order and number of moments (fixed)
N = 4;
Nmom = 35;
Nmom5 = 21;

% maximum number of time steps
nnmax = 20000000;
%nnmax = 5;

% initial correlation coefficients for joint Gaussian
r110 = 0.;
r101 = 0.;
r011 = 0.;

% largest dt to resolve collisions
dtmax = Kn;

%% shock problem %%%%%%%%%%%%
% initial densities
rhol = 1;
rhor = 0.01;

% initial mean velocities
U0 = 0;
V0 = 0;
W0 = 0;

% dimensionless temperature: used for scaling velocity so T=1 (do not change)
T = 1;

% set initial conditions to joint Gaussian with covariance
C200 = T;
C020 = T;
C002 = T;
C110 = r110*sqrt(C200*C020);
C101 = r101*sqrt(C200*C002);
C011 = r011*sqrt(C020*C002);

% initialize moments on "left" and "right"
Ml = InitializeM4_35(rhol,U0,V0,W0,C200,C110,C101,C020,C011,C002);
Mr = InitializeM4_35(rhor,U0,V0,W0,C200,C110,C101,C020,C011,C002);

%%
C200c = T;
C020c = T;
C002c = T;
C110c = r110*sqrt(C200*C020);
C101c = r101*sqrt(C200*C002);
C011c = r011*sqrt(C020*C002);
% magnitude of 3-D velocity = Ma
Uc = Ma/sqrt(2);
% initialize moments on "top" and "bottom" for crossing
Mt = InitializeM4_35(rhol,-Uc,-Uc,W0,C200c,C110c,C101c,C020c,C011c,C002c);
Mb = InitializeM4_35(rhol, Uc, Uc,W0,C200c,C110c,C101c,C020c,C011c,C002c);
%%

M = zeros(Np,Np,Nmom);
Mnp = M;
Mnpx = M;
Mnpy = M;
S = zeros(Np,Np,Nmom);
C = zeros(Np,Np,Nmom);
M5 = zeros(Np,Np,Nmom5);
S5 = zeros(Np,Np,Nmom5);
C5 = zeros(Np,Np,Nmom5);

%% initialize 35 3-D moments on 2-D spatial domain
% low-pressure background
for i = 1:Np
    for j = 1:Np
        for kk = 1:Nmom
            M(i,j,kk) = Mr(kk);
        end
    end
end
% high-pressure center (Csize = size of center region)
Csize = floor(0.1*Np) ;
Mint = Np/2 + 1;
Maxt = Np/2 + 1 + Csize;
Minb = Np/2 - Csize;
Maxb = Np/2;
for i = Minb:Maxb
    for j = Minb:Maxb
        for kk = 1:Nmom
            M(i,j,kk) = Mb(kk);
        end
    end
end
for i = Mint:Maxt
    for j = Mint:Maxt
        for kk = 1:Nmom
            M(i,j,kk) = Mt(kk);
        end
    end
end
%
% compute moments and plot initial conditions
for i = 1:Np
    for j = 1:Np
        MOM = zeros(Nmom,1);
        for kk=1:Nmom
            MOM(kk,1) = M(i,j,kk);
        end
        [CC,SS] = M2CS4_35(MOM);
        for kk=1:Nmom
            C(i,j,kk) = CC(kk);
            S(i,j,kk) = SS(kk);
        end
    end
end

for i = 1:Np
    for j = 1:Np
        MOM = zeros(Nmom,1);
        for kk=1:Nmom
            MOM(kk,1) = M(i,j,kk);
        end
        [MM5,CC5,SS5] = Moments5_3D(MOM);
        for kk=1:Nmom5
            M5(i,j,kk) = MM5(kk);
            C5(i,j,kk) = CC5(kk);
            S5(i,j,kk) = SS5(kk);
        end
    end
end

nmin = 1;
nmax = Np;

cc = 'k';

figure(1)
contour3D_plots
colormap sky

figure(2)
plot3Dsym_mom

figure(3)
plot3Dsym_C

figure(4)
plot3Dsym_S

pause(4)
%%

% name saved file
txt = ['riemann_3D_hyqmom35_crossing','_Np',num2str(Np),'_Kn',num2str(Kn),'_Ma',num2str(Ma),'.mat'];

%% time evolution begins here
t = 0.;
Fx = zeros(Np,Np,Nmom);
Fy = zeros(Np,Np,Nmom);
Mx = zeros(1,Nmom);  % closures for x flux
My = zeros(1,Nmom);  % closures for y flux
Mr = zeros(1,Nmom);  % realizable moments
vpxmin = zeros(Np,Np,1);
vpxmax = zeros(Np,Np,1);
vpymin = zeros(Np,Np,1);
vpymax = zeros(Np,Np,1);
v5xmin = zeros(Np,Np,1);
v5xmax = zeros(Np,Np,1);
v5ymin = zeros(Np,Np,1);
v5ymax = zeros(Np,Np,1);
v6xmin = zeros(Np,Np,1);
v6xmax = zeros(Np,Np,1);
v6ymin = zeros(Np,Np,1);
v6ymax = zeros(Np,Np,1);
nn = 0;

tic
while t<tmax && nn<nnmax
    nn = nn+1;
    
    % spatial fluxes, realizability checks, eigenvalues
    Mnp = M;
    parfor i = 1:Np
        for j = 1:Np
            MOM = zeros(Nmom,1);
            for kk = 1:Nmom
                MOM(kk) = M(i,j,kk) ;
            end
            % eigenvalues with hyperbolicity
            [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM,flag2D,Ma);
            [v6xmin(i,j),v6xmax(i,j),Mr] = eigenvalues6x_hyperbolic_3D(Mr,flag2D,Ma);
            [v6ymin(i,j),v6ymax(i,j),Mr] = eigenvalues6y_hyperbolic_3D(Mr,flag2D,Ma);
            [Mx,My,~,Mr] = Flux_closure35_and_realizable_3D(Mr,flag2D,Ma);
            % fluxes in the x direction
            Fx(i,j,:) = Mx;
            % fluxes in the y direction
            Fy(i,j,:) = My;
            % realizable moments
            Mnp(i,j,:)= Mr;
            %
            % compute eigenvalues for HLL
            % 1-D hyqmom for m500 eigenvalues in x direction
            MOM5 = [Mr(1) Mr(2) Mr(3) Mr(4) Mr(5)]; % m000 m100 m200 m300 m400
            [~,v5xmin(i,j),v5xmax(i,j)] = closure_and_eigenvalues(MOM5);
            %
            vpxmin(i,j)=min(v5xmin(i,j),v6xmin(i,j));
            vpxmax(i,j)=max(v5xmax(i,j),v6xmax(i,j));
            % 1-D hyqmom for m050 eigenvalues in y direction
            MOM5 = [Mr(1) Mr(6) Mr(10) Mr(13) Mr(15)]; % m000 m010 m020 m030 m040
            [~,v5ymin(i,j),v5ymax(i,j)] = closure_and_eigenvalues(MOM5);
            %
            vpymin(i,j)=min(v5ymin(i,j),v6ymin(i,j));
            vpymax(i,j)=max(v5ymax(i,j),v6ymax(i,j));
            % 
        end
    end
    M = Mnp;

    % fix time step based on largest eigenvalues in computational domain
    dt = CFL*dx/max([abs(vpxmax);abs(vpxmin);abs(vpymax);abs(vpymin)],[],'all');
    if t+dt>tmax
        dt=tmax-t;
    end
    dt = min(dt,dtmax);
    t = t+dt
    
    %% Euler for flux starts here
    % update moments due to spatial fluxes using method of lines and HLL
    Mnpx = M;
    parfor j = 1:Np
        VxMIN = zeros(Np,1);
        VxMAX = zeros(Np,1);
        MOM = zeros(Np,Nmom);
        FX = zeros(Np,Nmom);
        for i = 1:Np
            for kk = 1:Nmom
                MOM(i,kk)= M(i,j,kk);
                FX(i,kk) = Fx(i,j,kk);
            end
            VxMIN(i,1) = vpxmin(i,j);
            VxMAX(i,1) = vpxmax(i,j);
        end
        MNP = pas_HLL(MOM,FX,dt,dx,VxMIN,VxMAX);
        for i = 1:Np
            for kk = 1:Nmom
                Mnpx(i,j,kk) = MNP(i,kk);
            end
        end
    end
    %
    Mnpy = M;
    parfor i = 1:Np
        VyMIN = zeros(Np,1);
        VyMAX = zeros(Np,1);
        MOM = zeros(Np,Nmom);
        FY = zeros(Np,Nmom);
        for j = 1:Np
            for kk = 1:Nmom
                MOM(j,kk)= M(i,j,kk);
                FY(j,kk) = Fy(i,j,kk);
            end
            VyMIN(j,1) = vpymin(i,j);
            VyMAX(j,1) = vpymax(i,j);
        end
        MNP = pas_HLL(MOM,FY,dt,dy,VyMIN,VyMAX);
        for j = 1:Np
            for kk = 1:Nmom
                Mnpy(i,j,kk) = MNP(j,kk);
            end
        end
    end
    % end of Euler (NB: Mnp can have unrealizable moments)
    Mnp = Mnpx + Mnpy - M;
    %%
    M = Mnp;
    %
    % enforce realizability and hyperbolicity
    parfor i = 1:Np
        for j = 1:Np
            MOM = zeros(Nmom,1);
            for kk = 1:Nmom
                MOM(kk) = M(i,j,kk) ;
            end
            [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM,flag2D,Ma);
            [v6xmin(i,j),v6xmax(i,j),Mr] = eigenvalues6x_hyperbolic_3D(Mr,flag2D,Ma);
            [v6ymin(i,j),v6ymax(i,j),Mr] = eigenvalues6y_hyperbolic_3D(Mr,flag2D,Ma);
            [~,~,~,Mr] = Flux_closure35_and_realizable_3D(Mr,flag2D,Ma);
            % realizable moments
            Mnp(i,j,:)= Mr;
        end
    end
    M = Mnp;
    %
    % collision step using BGK
    parfor i = 1:Np
        for j = 1:Np
            MM = zeros(Nmom,1);
            for kk = 1:Nmom
                MM(kk,1) = M(i,j,kk);
            end
            MMC = collision35(MM,dt,Kn);
            for kk = 1:Nmom
                Mnp(i,j,kk) = MMC(kk,1);
            end
        end
    end
    M = Mnp;
    %
    % compute central and standardized moments, check 1-D realizability
    parfor i = 1:Np
        for j = 1:Np
            MOM = zeros(Nmom,1);
            for kk=1:Nmom
                MOM(kk,1) = M(i,j,kk);
            end
            [CC,SS] = M2CS4_35(MOM);
            for kk=1:Nmom
                C(i,j,kk) = CC(kk);
                S(i,j,kk) = SS(kk);
            end
        end
    end
    %
    if any(C(:,:,3) < 0,'all') 
        disp('pb C200 realizabilite apres pas temps')
        break
    end
    if any(C(:,:,10) < 0,'all')
        disp('pb C020 realizabilite apres pas temps')
        break
    end
    if any(C(:,:,20) < 0,'all')
        disp('pb C002 realizabilite apres pas temps')
        break
    end
    %
    if any(S(:,:,5)-1-S(:,:,4).^2 < 0,'all') 
        disp('pb H200 realizabilite apres pas temps')
        %break
    end
    if any(S(:,:,15)-1-S(:,:,13).^2 < 0,'all') 
        disp('pb H020 realizabilite apres pas temps')
        %break
    end
    if any(S(:,:,25)-1-S(:,:,23).^2 < 0,'all') 
        disp('pb H002 realizabilite apres pas temps')
        %break
    end

    figure(10)
    Cmoment3D_plots
    colormap sky
    pause(2)
    hold off

    [Diff,MaxDiff] = test_symmetry_2D(M,Np);
    MaxDiff 
end
toc

%% postprocessing for plots
for i = 1:Np
    for j = 1:Np
        MOM = zeros(Nmom,1);
        for kk=1:Nmom
            MOM(kk,1) = M(i,j,kk);
        end
        [MM5,CC5,SS5] = Moments5_3D(MOM);
        for kk=1:Nmom5
            M5(i,j,kk) = MM5(kk);
            C5(i,j,kk) = CC5(kk);
            S5(i,j,kk) = SS5(kk);
        end
    end
end

lam6xa = zeros(Np,Np,6);
lam6xb = zeros(Np,Np,6);
lam6ya = zeros(Np,Np,6);
lam6yb = zeros(Np,Np,6);
for i = 1:Np
    for j = 1:Np
        M1 = zeros(Nmom,1);
        for kk=1:Nmom
            M1(kk,1) = M(i,j,kk);
        end
        %
        m000 = M1(1);
        m100 = M1(2);
        m200 = M1(3);
        m300 = M1(4);
        m400 = M1(5);
        m010 = M1(6);
        m110 = M1(7);
        m210 = M1(8);
        m310 = M1(9);
        m020 = M1(10);
        m120 = M1(11);
        m220 = M1(12);
        m030 = M1(13);
        m130 = M1(14);
        m040 = M1(15);
        m001 = M1(16);
        m101 = M1(17);
        m201 = M1(18);
        m301 = M1(19);
        m002 = M1(20);
        m102 = M1(21);
        m202 = M1(22);
        m003 = M1(23);
        m103 = M1(24);
        m004 = M1(25);
        m011 = M1(26);
        m021 = M1(29);
        m031 = M1(31);
        m012 = M1(32);
        m013 = M1(34);
        m022 = M1(35);
        J6 = jacobian6(m000,m010,m020,m030,m040,m100,m110,m120,m130,m200,m210,m220,m300,m310,m400);
        lam6xa(i,j,:) = eig(J6);
        J6 = jacobian6(m000,m001,m002,m003,m004,m100,m101,m102,m103,m200,m201,m202,m300,m301,m400);
        lam6xb(i,j,:) = eig(J6);
        J6 = jacobian6(m000,m100,m200,m300,m400,m010,m110,m210,m310,m020,m120,m220,m030,m130,m040);
        lam6ya(i,j,:) = eig(J6);
        J6 = jacobian6(m000,m001,m002,m003,m004,m010,m011,m012,m013,m020,m021,m022,m030,m031,m040);
        lam6yb(i,j,:) = eig(J6);
    end
end

% save date
save(txt)

% Also save for comparison with new implementation
fprintf('\n=== Original Implementation Results ===\n');
fprintf('Grid: %dx%d\n', Np, Np);
fprintf('Final time: %.6f\n', t);
fprintf('Number of timesteps: %d\n', nn);
fprintf('M(20,20,1) = %.15e\n', M(20,20,1));
fprintf('M(20,20,2) = %.15e\n', M(20,20,2));

original_results = struct();
original_results.M = M;
original_results.C = C;
original_results.S = S;
original_results.grid.xm = xm;
original_results.grid.ym = ym;
original_results.params.Np = Np;
original_results.params.tmax = tmax;
original_results.params.final_time = t;
original_results.params.num_timesteps = nn;
save('test_original_Np40_results.mat', 'original_results');
fprintf('Saved to test_original_Np40_results.mat\n');

return;  % Skip plots for automated testing

%% plots
nmin = 1;
nmax = Np;

cc = 'r';

figure(2)
%nc = 5; nl = 4;
plot3Dsym_mom

figure(3)
plot3Dsym_C

figure(4)
plot3Dsym_S

figure(5)
Y1 = 0*xm;
Y2 = 0*xm;
Y3 = 0*xm;
Y4 = 0*xm;
for i=1:Np
    Y1(i) = v5xmin(i,i);
    Y2(i) = v6xmin(i,i);
    Y3(i) = v5xmax(i,i);
    Y4(i) = v6xmax(i,i);
end
plot(xm,Y1,'k',xm,Y2,'r',xm,Y3,'k',xm,Y4,'r')

figure(6)
LAMXa = zeros(Np,6);
LAMXb = zeros(Np,6);
for i = 1:Np
    for kk = 1:6
        LAMXa(i,kk) = real(lam6xa(i,i,kk));
        LAMXb(i,kk) = real(lam6xb(i,i,kk));
    end
end
plot(xm,LAMXa,'o',xm,LAMXb,'p')

figure(7)
Y1 = 0*xm;
Y2 = 0*xm;
Y3 = 0*xm;
Y4 = 0*xm;
for i=1:Np
    Y1(i) = v5ymin(i,i);
    Y2(i) = v6ymin(i,i);
    Y3(i) = v5ymax(i,i);
    Y4(i) = v6ymax(i,i);
end
plot(xm,Y1,'k',xm,Y2,'r',xm,Y3,'k',xm,Y4,'r')

figure(8)
LAMYa = zeros(Np,6);
LAMYb = zeros(Np,6);
for j = 1:Np
    for kk = 1:6
        LAMYa(j,kk) = real(lam6ya(j,j,kk));
        LAMYb(j,kk) = real(lam6yb(j,j,kk));
    end
end
plot(ym,LAMYa,'o',ym,LAMYb,'p')

figure(9)
contour3D_plots
colormap sky

figure(10)
Cmoment3D_plots
colormap sky

figure(11)
Smoment3D_plots
colormap sky

figure(12)
hyperbolic3D_plots
colormap sky