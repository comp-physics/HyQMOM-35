"""
Visualization module for HyQMOM simulation results.

This module provides PyPlot-based plotting functions equivalent to MATLAB's
simulation_plots.m, for visualizing final simulation results.
"""

using PyPlot
using ColorSchemes
using Printf

# Import required functions from HyQMOM modules
using LinearAlgebra: det

"""
    plot_final_results(M, xm, ym, Np, Nmom; save_figures=false, output_dir=".")

Create comprehensive visualization of final simulation results.

This is the main plotting function that creates all final result figures
(Figures 2-12 in MATLAB code).

# Arguments
- `M`: Moment array (Np x Np x Nmom)
- `xm`: X-coordinates of cell centers
- `ym`: Y-coordinates of cell centers
- `Np`: Grid size
- `Nmom`: Number of moments (35)
- `save_figures`: Save figures to disk (default: false)
- `output_dir`: Directory for saved figures (default: ".")

# Figures Created
- Figure 2: Moment line plots along diagonal
- Figure 3: Central moment line plots
- Figure 4: Standardized moment line plots  
- Figure 9: Contour plots (12 panels)
- Figure 10: C-moment contour plots (16 panels)
- Figure 11: S-moment contour plots (12 panels)
- Figure 12: Hyperbolicity/eigenvalue plots (9 panels)
"""
function plot_final_results(M, xm, ym, Np, Nmom; save_figures=false, output_dir=".")
    println("\n" * "="^70)
    println("Generating visualization of final results...")
    println("="^70)
    
    # Compute C and S moments from M
    println("Computing central and standardized moments...")
    C = zeros(Np, Np, Nmom)
    S = zeros(Np, Np, Nmom)
    for i = 1:Np
        for j = 1:Np
            MOM = M[i, j, :]
            CC, SS = M2CS4_35(MOM)
            C[i, j, :] = CC
            S[i, j, :] = SS
        end
    end
    
    # Compute 5th order moments
    println("Computing 5th order moments...")
    Nmom5 = 3  # C500, C410, C320
    M5 = zeros(Np, Np, Nmom5)
    C5 = zeros(Np, Np, Nmom5)
    S5 = zeros(Np, Np, Nmom5)
    for i = 1:Np
        for j = 1:Np
            MOM = M[i, j, :]
            MM5, CC5, SS5 = Moments5_3D(MOM)
            M5[i, j, :] = MM5
            C5[i, j, :] = CC5
            S5[i, j, :] = SS5
        end
    end
    
    # Compute eigenvalues for hyperbolicity plots
    println("Computing eigenvalues...")
    eig_data = compute_grid_eigenvalues(M, Np, Nmom)
    
    # Create plots
    nmin = 1
    nmax = Np
    cc = "r"  # Red color for final results
    
    println("Creating Figure 2: Moment line plots...")
    figure(2)
    clf()
    plot_3Dsym_moments(xm, M, nmin, nmax, cc, Np)
    tight_layout()
    if save_figures
        savefig(joinpath(output_dir, "fig02_moments.png"), dpi=150)
    end
    
    println("Creating Figure 3: Central moment line plots...")
    figure(3)
    clf()
    plot_3Dsym_central(xm, M, C, C5, nmin, nmax, cc, Np)
    tight_layout()
    if save_figures
        savefig(joinpath(output_dir, "fig03_central_moments.png"), dpi=150)
    end
    
    println("Creating Figure 4: Standardized moment line plots...")
    figure(4)
    clf()
    plot_3Dsym_standardized(xm, S, S5, nmin, nmax, cc, Np)
    tight_layout()
    if save_figures
        savefig(joinpath(output_dir, "fig04_standardized_moments.png"), dpi=150)
    end
    
    println("Creating Figure 9: Contour plots...")
    figure(9)
    clf()
    contour_plots_3D(xm, ym, M, C, S, Np)
    tight_layout()
    if save_figures
        savefig(joinpath(output_dir, "fig09_contours.png"), dpi=150)
    end
    
    println("Creating Figure 10: C-moment contour plots...")
    figure(10)
    clf()
    Cmoment_plots_3D(xm, ym, S, Np)
    tight_layout()
    if save_figures
        savefig(joinpath(output_dir, "fig10_cmoments.png"), dpi=150)
    end
    
    println("Creating Figure 11: S-moment contour plots...")
    figure(11)
    clf()
    Smoment_plots_3D(xm, ym, S, Np)
    tight_layout()
    if save_figures
        savefig(joinpath(output_dir, "fig11_smoments.png"), dpi=150)
    end
    
    println("Creating Figure 12: Hyperbolicity plots...")
    figure(12)
    clf()
    hyperbolic_plots_3D(xm, ym, eig_data, Np)
    tight_layout()
    if save_figures
        savefig(joinpath(output_dir, "fig12_hyperbolicity.png"), dpi=150)
    end
    
    println("="^70)
    println("Visualization complete! Created 7 figures.")
    if save_figures
        println("Figures saved to: $output_dir")
        println("  - fig02_moments.png")
        println("  - fig03_central_moments.png")
        println("  - fig04_standardized_moments.png")
        println("  - fig09_contours.png")
        println("  - fig10_cmoments.png")
        println("  - fig11_smoments.png")
        println("  - fig12_hyperbolicity.png")
    else
        println("Figures displayed in windows.")
        println("Call PyPlot.show() to display if not visible, or close windows to continue.")
        # Show the figures
        PyPlot.show()
    end
    println("="^70)
    
    return nothing
end

"""
    contour_plots_3D(xm, ym, M, C, S, Np)

Create 12-panel contour plot figure showing density, velocities, moments,
and realizability quantities.

Equivalent to MATLAB's contour_plots_3D function (simulation_plots.m lines 173-294).
"""
function contour_plots_3D(xm, ym, M, C, S, Np)
    # Set colormap to match MATLAB 'sky'
    set_sky_colormap()
    
    X, Y = meshgrid(xm, ym)
    
    # Panel 1: Density
    subplot(3, 4, 1)
    Z = M[:, :, 1]
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("Density")
    colorbar()
    
    # Panel 2: U velocity
    subplot(3, 4, 2)
    Z = M[:, :, 2] ./ M[:, :, 1]
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("U velocity")
    colorbar()
    
    # Panel 3: V velocity
    subplot(3, 4, 3)
    Z = M[:, :, 6] ./ M[:, :, 1]
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("V velocity")
    colorbar()
    
    # Panel 4: S_110
    subplot(3, 4, 4)
    Z = S[:, :, 7]
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("S_{110}")
    colorbar()
    
    # Panel 5: C_200
    subplot(3, 4, 5)
    Z = C[:, :, 3]
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("C_{200}")
    colorbar()
    
    # Panel 6: C_020
    subplot(3, 4, 6)
    Z = C[:, :, 10]
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("C_{020}")
    colorbar()
    
    # Panel 7: C_002
    subplot(3, 4, 7)
    Z = C[:, :, 20]
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("C_{002}")
    colorbar()
    
    # Panel 8: |Delta_1|
    subplot(3, 4, 8)
    Z = 1 .+ 2 .* S[:, :, 7] .* S[:, :, 17] .* S[:, :, 26] .- 
        S[:, :, 7].^2 .- S[:, :, 17].^2 .- S[:, :, 26].^2
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("|Δ_1|")
    colorbar()
    
    # Panel 9: H_200
    subplot(3, 4, 9)
    Z = max.(0, S[:, :, 5] .- S[:, :, 4].^2 .- 1)
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("H_{200}")
    colorbar()
    
    # Panel 10: H_020
    subplot(3, 4, 10)
    Z = max.(0, S[:, :, 15] .- S[:, :, 13].^2 .- 1)
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("H_{020}")
    colorbar()
    
    # Panel 11: H_002
    subplot(3, 4, 11)
    Z = max.(0, S[:, :, 25] .- S[:, :, 23].^2 .- 1)
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("H_{002}")
    colorbar()
    
    # Panel 12: |Delta_2*|
    subplot(3, 4, 12)
    dDel2 = zeros(Np, Np)
    for i = 1:Np
        for j = 1:Np
            # Extract standardized moments
            S300 = S[i, j, 4]; S400 = S[i, j, 5]
            S110 = S[i, j, 7]; S210 = S[i, j, 8]; S310 = S[i, j, 9]
            S120 = S[i, j, 11]; S220 = S[i, j, 12]
            S030 = S[i, j, 13]; S130 = S[i, j, 14]; S040 = S[i, j, 15]
            S101 = S[i, j, 17]; S201 = S[i, j, 18]; S301 = S[i, j, 19]
            S102 = S[i, j, 21]; S202 = S[i, j, 22]
            S003 = S[i, j, 23]; S103 = S[i, j, 24]; S004 = S[i, j, 25]
            S011 = S[i, j, 26]; S111 = S[i, j, 27]
            S211 = S[i, j, 28]; S021 = S[i, j, 29]; S121 = S[i, j, 30]
            S031 = S[i, j, 31]; S012 = S[i, j, 32]
            S112 = S[i, j, 33]; S013 = S[i, j, 34]; S022 = S[i, j, 35]
            
            D2 = delta2star3D(S300, S400, S110, S210, S310, S120, S220, S030, S130, S040,
                              S101, S201, S301, S102, S202, S003, S103, S004,
                              S011, S111, S211, S021, S121, S031, S012, S112, S013, S022)
            
            # Compute minimum of principal minors
            dD2s = minimum([det(D2[1:1, 1:1]), det(D2[1:2, 1:2]), det(D2[1:3, 1:3]),
                            det(D2[1:4, 1:4]), det(D2[1:5, 1:5]), det(D2)])
            dDel2[i, j] = max(0, dD2s)
        end
    end
    Z = dDel2
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("|Δ_2*|")
    colorbar()
end

"""
    Cmoment_plots_3D(xm, ym, S, Np)

Create 16-panel C-moment contour plot figure.

Equivalent to MATLAB's Cmoment_plots_3D function (simulation_plots.m lines 296-426).
"""
function Cmoment_plots_3D(xm, ym, S, Np)
    set_sky_colormap()
    X, Y = meshgrid(xm, ym)
    
    # Define the 16 moments to plot with their indices and labels
    moments = [
        (17, "S_{101}"), (26, "S_{011}"), (18, "S_{201}"), (29, "S_{021}"),
        (21, "S_{102}"), (32, "S_{012}"), (22, "S_{202}"), (35, "S_{022}"),
        (24, "S_{103}"), (34, "S_{013}"), (19, "S_{301}"), (31, "S_{031}"),
        (27, "S_{111}"), (28, "S_{211}"), (30, "S_{121}"), (33, "S_{112}")
    ]
    
    for (panel_idx, (moment_idx, label)) in enumerate(moments)
        subplot(4, 4, panel_idx)
        
        if moment_idx == 22  # S_202 - 1
            Z = S[:, :, moment_idx] .- 1
        elseif moment_idx == 35  # S_022 - 1
            Z = S[:, :, moment_idx] .- 1
        elseif moment_idx == 33  # S_112 - S_110
            Z = real.(S[:, :, moment_idx] .- S[:, :, 7])
        else
            Z = S[:, :, moment_idx]
        end
        
        contourf(X, Y, Z', 50, cmap="sky")
        axis("square")
        
        # Adjust title for special cases
        if moment_idx == 22
            title("S_{202}-1")
        elseif moment_idx == 35
            title("S_{022}-1")
        elseif moment_idx == 33
            title("S_{112}-S_{110}")
        else
            title(label)
        end
        colorbar()
    end
end

"""
    Smoment_plots_3D(xm, ym, S, Np)

Create 12-panel S-moment contour plot figure.

Equivalent to MATLAB's Smoment_plots_3D function (simulation_plots.m lines 428-528).
"""
function Smoment_plots_3D(xm, ym, S, Np)
    set_sky_colormap()
    X, Y = meshgrid(xm, ym)
    
    cmin = 0
    cmax = 100
    
    # Define the 12 moments to plot
    moments = [
        (4, "S_{300}"), (8, "S_{210}"), (11, "S_{120}"), (32, "S_{012}"),
        (13, "S_{030}"), (5, "S_{400}"), (9, "S_{310}"), (12, "S_{220}"),
        (14, "S_{130}"), (25, "S_{004}"), (33, "S_{112}"), (35, "S_{022}")
    ]
    
    for (panel_idx, (moment_idx, label)) in enumerate(moments)
        subplot(3, 4, panel_idx)
        
        Z = S[:, :, moment_idx]
        
        # Apply clipping for first panel (S_300)
        if moment_idx == 4
            Z = max.(-cmax, min.(cmax, Z))
        end
        
        contourf(X, Y, Z', 50, cmap="sky")
        axis("square")
        title(label)
        colorbar()
    end
end

"""
    hyperbolic_plots_3D(xm, ym, eig_data, Np)

Create 9-panel eigenvalue contour plot figure for hyperbolicity analysis.

Equivalent to MATLAB's hyperbolic_plots_3D function (simulation_plots.m lines 530-574).
"""
function hyperbolic_plots_3D(xm, ym, eig_data, Np)
    set_sky_colormap()
    X, Y = meshgrid(xm, ym)
    
    # Extract eigenvalues for UV plane (X-direction)
    lam6x = eig_data[:lam6x]
    
    # Plot first 6 eigenvalues
    for k = 1:6
        subplot(3, 3, k)
        Z = lam6x[:, :, k]
        contourf(X, Y, Z', 50, cmap="sky")
        axis("square")
        title("λ_$k")
        colorbar()
    end
    
    # Plot eigenvalue differences
    subplot(3, 3, 7)
    Z = lam6x[:, :, 2] .- lam6x[:, :, 1]
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("λ_2 - λ_1")
    colorbar()
    
    subplot(3, 3, 8)
    Z = lam6x[:, :, 4] .- lam6x[:, :, 3]
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("λ_4 - λ_3")
    colorbar()
    
    subplot(3, 3, 9)
    Z = lam6x[:, :, 6] .- lam6x[:, :, 5]
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("λ_6 - λ_5")
    colorbar()
end

"""
    plot_3Dsym_moments(xm, M, nmin, nmax, cc, Np)

Create 12-panel line plot of moments vs x along diagonal (y=x).

Equivalent to MATLAB's plot_3Dsym_moments function (simulation_plots.m lines 576-700).
"""
function plot_3Dsym_moments(xm, M, nmin, nmax, cc, Np)
    nc = 4
    nl = 3
    
    # Define moments to plot: (moment_idx, label)
    moments = [
        (1, "M_{000}"), (2, "M_{100}"), (3, "M_{200}"), (20, "M_{002}"),
        (4, "M_{300}"), (5, "M_{400}"), (7, "M_{110}"), (35, "M_{022}"),
        (8, "M_{210}"), (9, "M_{310}"), (12, "M_{220}"), (25, "M_{004}")
    ]
    
    for (panel_idx, (moment_idx, label)) in enumerate(moments)
        subplot(nl, nc, panel_idx)
        hold(true)
        
        # Extract diagonal values
        Y = [M[i, i, moment_idx] for i = 1:Np]
        
        plot(xm[nmin:nmax], Y[nmin:nmax], color=cc, linewidth=2)
        xlabel("x")
        ylabel(label)
        grid(true)
        ax = gca()
        ax.tick_params(labelsize=18)
    end
end

"""
    plot_3Dsym_central(xm, M, C, C5, nmin, nmax, cc, Np)

Create 12-panel line plot of central moments vs x along diagonal.

Equivalent to MATLAB's plot_3Dsym_central function (simulation_plots.m lines 702-827).
"""
function plot_3Dsym_central(xm, M, C, C5, nmin, nmax, cc, Np)
    nl = 3
    nc = 4
    
    # Panel 1: M_000 (density)
    subplot(nl, nc, 1)
    hold(true)
    Y = [M[i, i, 1] for i = 1:Np]
    plot(xm[nmin:nmax], Y[nmin:nmax], color=cc, linewidth=2)
    xlabel("x")
    ylabel("M_{000}")
    grid(true)
    ax = gca()
    ax.tick_params(labelsize=18)
    
    # Panel 2: u_1 (mean velocity)
    subplot(nl, nc, 2)
    hold(true)
    umean = M[:, :, 2] ./ M[:, :, 1]
    Y = [umean[i, i] for i = 1:Np]
    plot(xm[nmin:nmax], Y[nmin:nmax], color=cc, linewidth=2)
    xlabel("x")
    ylabel("u_1")
    grid(true)
    ax = gca()
    ax.tick_params(labelsize=18)
    
    # Central moments C
    c_moments = [
        (3, "C_{200}"), (4, "C_{300}"), (5, "C_{400}"), (7, "C_{110}"),
        (8, "C_{210}"), (9, "C_{310}"), (12, "C_{220}")
    ]
    
    for (panel_idx, (c_idx, label)) in enumerate(c_moments)
        subplot(nl, nc, panel_idx + 2)
        hold(true)
        Y = [C[i, i, c_idx] for i = 1:Np]
        plot(xm[nmin:nmax], Y[nmin:nmax], color=cc, linewidth=2)
        xlabel("x")
        ylabel(label)
        grid(true)
        ax = gca()
        ax.tick_params(labelsize=18)
    end
    
    # 5th order central moments
    c5_moments = [(1, "C_{500}"), (2, "C_{410}"), (3, "C_{320}")]
    
    for (panel_idx, (c5_idx, label)) in enumerate(c5_moments)
        subplot(nl, nc, panel_idx + 9)
        hold(true)
        Y = [C5[i, i, c5_idx] for i = 1:Np]
        plot(xm[nmin:nmax], Y[nmin:nmax], color=cc, linewidth=2)
        xlabel("x")
        ylabel(label)
        grid(true)
        ax = gca()
        ax.tick_params(labelsize=18)
    end
end

"""
    plot_3Dsym_standardized(xm, S, S5, nmin, nmax, cc, Np)

Create 12-panel line plot of standardized moments and realizability quantities.

Equivalent to MATLAB's plot_3Dsym_standardized function (simulation_plots.m lines 829-964).
"""
function plot_3Dsym_standardized(xm, S, S5, nmin, nmax, cc, Np)
    nl = 3
    nc = 4
    
    # Panel 1: S_300
    subplot(nl, nc, 1)
    hold(true)
    Y = [S[i, i, 4] for i = 1:Np]
    plot(xm[nmin:nmax], Y[nmin:nmax], color=cc, linewidth=2)
    xlabel("x")
    ylabel("S_{300}")
    grid(true)
    ax = gca()
    ax.tick_params(labelsize=18)
    
    # Panel 2: S_400
    subplot(nl, nc, 2)
    hold(true)
    Y = [S[i, i, 5] for i = 1:Np]
    plot(xm[nmin:nmax], Y[nmin:nmax], color=cc, linewidth=2)
    xlabel("x")
    ylabel("S_{400}")
    grid(true)
    ax = gca()
    ax.tick_params(labelsize=18)
    
    # Panel 3: H_400
    subplot(nl, nc, 3)
    hold(true)
    H40 = S[:, :, 5] .- S[:, :, 4].^2 .- 1
    Y = [H40[i, i] for i = 1:Np]
    plot(xm[nmin:nmax], Y[nmin:nmax], color=cc, linewidth=2)
    xlabel("x")
    ylabel("H_{400}")
    grid(true)
    ax = gca()
    ax.tick_params(labelsize=18)
    
    # Panel 4: S_500
    subplot(nl, nc, 4)
    hold(true)
    Y = [S5[i, i, 1] for i = 1:Np]
    plot(xm[nmin:nmax], Y[nmin:nmax], color=cc, linewidth=2)
    xlabel("x")
    ylabel("S_{500}")
    grid(true)
    ax = gca()
    ax.tick_params(labelsize=18)
    
    # Panels 5-9: S moments
    s_moments = [
        (7, "S_{110}"), (8, "S_{210}"), (9, "S_{310}"), (12, "S_{220}")
    ]
    
    for (panel_idx, (s_idx, label)) in enumerate(s_moments)
        subplot(nl, nc, panel_idx + 4)
        hold(true)
        Y = [S[i, i, s_idx] for i = 1:Np]
        plot(xm[nmin:nmax], Y[nmin:nmax], color=cc, linewidth=2)
        xlabel("x")
        ylabel(label)
        grid(true)
        ax = gca()
        ax.tick_params(labelsize=18)
    end
    
    # Panel 9: S_410
    subplot(nl, nc, 9)
    hold(true)
    Y = [S5[i, i, 2] for i = 1:Np]
    plot(xm[nmin:nmax], Y[nmin:nmax], color=cc, linewidth=2)
    xlabel("x")
    ylabel("S_{410}")
    grid(true)
    ax = gca()
    ax.tick_params(labelsize=18)
    
    # Panel 10: S_320
    subplot(nl, nc, 10)
    hold(true)
    Y = [S5[i, i, 3] for i = 1:Np]
    plot(xm[nmin:nmax], Y[nmin:nmax], color=cc, linewidth=2)
    xlabel("x")
    ylabel("S_{320}")
    grid(true)
    ax = gca()
    ax.tick_params(labelsize=18)
    
    # Panel 11: |Delta_1|
    subplot(nl, nc, 11)
    hold(true)
    Y = [1 + 2*S[i, i, 7]*S[i, i, 17]*S[i, i, 26] - 
         S[i, i, 7]^2 - S[i, i, 17]^2 - S[i, i, 26]^2 for i = 1:Np]
    plot(xm[nmin:nmax], Y[nmin:nmax], color=cc, linewidth=2)
    xlabel("x")
    ylabel("|Δ_1|")
    grid(true)
    ax = gca()
    ax.tick_params(labelsize=18)
    
    # Panel 12: |Delta_2*|
    subplot(nl, nc, 12)
    hold(true)
    Y = zeros(Np)
    for i = 1:Np
        S300 = S[i, i, 4]; S400 = S[i, i, 5]; S110 = S[i, i, 7]
        S210 = S[i, i, 8]; S310 = S[i, i, 9]
        S120 = S[i, i, 11]; S220 = S[i, i, 12]; S030 = S[i, i, 13]
        S130 = S[i, i, 14]; S040 = S[i, i, 15]
        S101 = S[i, i, 17]; S201 = S[i, i, 18]; S301 = S[i, i, 19]
        S102 = S[i, i, 21]; S202 = S[i, i, 22]
        S003 = S[i, i, 23]; S103 = S[i, i, 24]; S004 = S[i, i, 25]
        S011 = S[i, i, 26]; S111 = S[i, i, 27]
        S211 = S[i, i, 28]; S021 = S[i, i, 29]; S121 = S[i, i, 30]
        S031 = S[i, i, 31]; S012 = S[i, i, 32]
        S112 = S[i, i, 33]; S013 = S[i, i, 34]; S022 = S[i, i, 35]
        
        D2 = delta2star3D(S300, S400, S110, S210, S310, S120, S220, S030, S130, S040,
                          S101, S201, S301, S102, S202, S003, S103, S004,
                          S011, S111, S211, S021, S121, S031, S012, S112, S013, S022)
        Y[i] = det(D2)
    end
    plot(xm[nmin:nmax], Y[nmin:nmax], color=cc, linewidth=2)
    xlabel("x")
    ylabel("|Δ_2*|")
    grid(true)
    ax = gca()
    ax.tick_params(labelsize=18)
end

# ============================================================================
# Helper Functions
# ============================================================================

"""
    meshgrid(x, y)

Create 2D coordinate matrices from coordinate vectors (MATLAB equivalent).
"""
function meshgrid(x, y)
    X = repeat(reshape(x, 1, :), length(y), 1)
    Y = repeat(y, 1, length(x))
    return X, Y
end

"""
    set_sky_colormap()

Set colormap to match MATLAB's 'sky' colormap.
"""
function set_sky_colormap()
    # Register custom 'sky' colormap if not already done
    # MATLAB 'sky' goes from dark blue to light blue/cyan
    try
        PyPlot.register_cmap("sky", ColorSchemes.ice)
    catch
        # If registration fails, just use a similar blue colormap
        nothing
    end
end

"""
    compute_grid_eigenvalues(M, Np, Nmom)

Compute eigenvalues for all grid points for hyperbolicity analysis.

Returns a dictionary with eigenvalue arrays for different directions.
"""
function compute_grid_eigenvalues(M, Np, Nmom)
    lam6x = zeros(Np, Np, 6)
    lam6y = zeros(Np, Np, 6)
    lam6z = zeros(Np, Np, 6)
    lam6w = zeros(Np, Np, 6)
    
    for i = 1:Np
        for j = 1:Np
            MOM = M[i, j, :]
            
            # X-direction: UV plane
            M_slice = axis_moment_slice(MOM, 1)
            J6 = jacobian6(M_slice...)
            eigs = eigvals(J6)
            lam6x[i, j, :] = sort(real.(eigs))
            
            # Y-direction: VU plane
            M_slice = axis_moment_slice(MOM, 2)
            J6 = jacobian6(M_slice...)
            eigs = eigvals(J6)
            lam6y[i, j, :] = sort(real.(eigs))
            
            # X-direction: UW plane
            M_slice = axis_moment_slice(MOM, 3)
            J6 = jacobian6(M_slice...)
            eigs = eigvals(J6)
            lam6z[i, j, :] = sort(real.(eigs))
            
            # Y-direction: VW plane
            M_slice = axis_moment_slice(MOM, 4)
            J6 = jacobian6(M_slice...)
            eigs = eigvals(J6)
            lam6w[i, j, :] = sort(real.(eigs))
        end
    end
    
    return Dict(
        :lam6x => lam6x,
        :lam6y => lam6y,
        :lam6z => lam6z,
        :lam6w => lam6w
    )
end

