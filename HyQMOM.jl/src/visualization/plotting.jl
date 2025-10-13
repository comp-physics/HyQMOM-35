"""
Visualization module for HyQMOM simulation results.

This module provides PyPlot-based plotting functions equivalent to MATLAB's
simulation_plots.m, for visualizing final simulation results.
"""

using PyPlot
using PyCall
using ColorSchemes
using Printf

# Import required functions from HyQMOM modules
using LinearAlgebra: det

"""
    format_colorbar(; shrink=0.8, aspect=20)

Create a properly sized colorbar with formatted tick labels.

# Arguments
- `shrink`: Scale factor for colorbar size (default: 0.8)
- `aspect`: Aspect ratio for colorbar (default: 20, makes it narrower)

# Returns
- Colorbar object with formatted labels
"""
function format_colorbar(; shrink=0.8, aspect=20)
    cbar = colorbar(shrink=shrink, aspect=aspect)
    # Format tick labels to use scientific notation with 2 decimal places
    cbar.formatter.set_powerlimits((-2, 2))
    cbar.formatter.set_useMathText(true)
    cbar.ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useMathText=true))
    cbar.ax.ticklabel_format(style="scientific", axis="y", scilimits=(-2, 2))
    # Reduce number of ticks for cleaner appearance
    cbar.ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=5))
    cbar.update_ticks()
    return cbar
end

"""
    plot_final_results(M, xm, ym, Np, Nmom; save_figures=false, output_dir=".", zm=nothing, Nz=1)

Create comprehensive visualization of final simulation results.

This is the main plotting function that creates all final result figures
(Figures 2-12 in MATLAB code).

# Arguments
- `M`: Moment array (Np x Np x Nmom) for 2D or (Np x Np x Nz x Nmom) for 3D
- `xm`: X-coordinates of cell centers
- `ym`: Y-coordinates of cell centers
- `Np`: Grid size
- `Nmom`: Number of moments (35)
- `save_figures`: Save figures to disk (default: false)
- `output_dir`: Directory for saved figures (default: ".")
- `zm`: Z-coordinates of cell centers (optional, for 3D)
- `Nz`: Number of z grid points (default: 1)

# Figures Created
- Figure 2: Moment line plots along diagonal
- Figure 3: Central moment line plots
- Figure 4: Standardized moment line plots  
- Figure 9: Contour plots (12 panels)
- Figure 10: C-moment contour plots (16 panels)
- Figure 11: S-moment contour plots (12 panels)
- Figure 12: Hyperbolicity/eigenvalue plots (9 panels)

# Notes
For 3D simulations (Nz > 1), plots the middle z-slice.
"""
function plot_final_results(M, xm, ym, Np, Nmom; save_figures=false, output_dir=".", zm=nothing, Nz=1)
    println("\n" * "="^70)
    println("Generating visualization of final results...")
    
    # Handle 3D arrays
    if ndims(M) == 4
        Nz = size(M, 3)
        
        # For true 3D simulations (Nz > 1), create comprehensive 3D visualizations
        if Nz > 1 && zm !== nothing
            println("3D simulation detected (Nz=$Nz): Creating comprehensive 3D visualizations...")
            plot_3d_diagnostics(M, xm, ym, zm, Np; Nz=Nz, Nmom=Nmom, save_figures=save_figures, output_dir=output_dir)
            
            # Also create standard 2D plots for middle z-slice
            k_mid = div(Nz, 2) + 1
            println("\nAlso creating standard 2D plots for middle z-slice $k_mid (z = $(zm[k_mid]))...")
            M_plot = M[:, :, k_mid, :]
        else
            # Nz=1 or no z-coordinates provided: treat as 2D
            k_mid = div(Nz, 2) + 1
            if zm !== nothing
                println("Extracting z-slice $k_mid of $Nz (z = $(zm[k_mid]))")
            else
                println("3D array with Nz=$Nz")
            end
            M_plot = M[:, :, k_mid, :]
        end
    else
        M_plot = M  # Already 2D
        println("2D simulation")
    end
    
    println("="^70)
    
    # Set up the custom sky colormap
    set_sky_colormap()
    
    # Check for NaNs/Infs in M_plot
    if any(isnan.(M_plot)) || any(isinf.(M_plot))
        @warn "M_plot contains NaN or Inf values - skipping 2D standard plots"
        return
    end
    
    # Compute C and S moments from M_plot
    println("Computing central and standardized moments...")
    C = zeros(Np, Np, Nmom)
    S = zeros(Np, Np, Nmom)
    for i = 1:Np
        for j = 1:Np
            MOM = M_plot[i, j, :]
            CC, SS = M2CS4_35(MOM)
            C[i, j, :] = CC
            S[i, j, :] = SS
        end
    end
    
    # Compute 5th order moments
    println("Computing 5th order moments...")
    Nmom5 = 21  # Full 5th order moments (21 elements)
    M5 = zeros(Np, Np, Nmom5)
    C5 = zeros(Np, Np, Nmom5)
    S5 = zeros(Np, Np, Nmom5)
    for i = 1:Np
        for j = 1:Np
            MOM = M_plot[i, j, :]
            MM5, CC5, SS5 = Moments5_3D(MOM)
            M5[i, j, :] = MM5
            C5[i, j, :] = CC5
            S5[i, j, :] = SS5
        end
    end
    
    # Compute eigenvalues for hyperbolicity plots
    println("Computing eigenvalues...")
    eig_data = compute_grid_eigenvalues(M_plot, Np, Nmom)
    
    # Create plots
    nmin = 1
    nmax = Np
    cc = "r"  # Red color for final results
    
    println("Creating Figure 1: Moment line plots...")
    fig1 = figure(1)
    clf()
    plot_3Dsym_moments(xm, M_plot, nmin, nmax, cc, Np)
    tight_layout()
    if save_figures
        savefig(joinpath(output_dir, "fig01_moments.png"), dpi=150)
    end
    fig1.canvas.draw()  # Force drawing to ensure figure persists
    
    println("Creating Figure 2: Central moment line plots...")
    fig2 = figure(2)
    clf()
    plot_3Dsym_central(xm, M_plot, C, C5, nmin, nmax, cc, Np)
    tight_layout()
    if save_figures
        savefig(joinpath(output_dir, "fig02_central_moments.png"), dpi=150)
    end
    fig2.canvas.draw()  # Force drawing to ensure figure persists
    
    println("Creating Figure 3: Standardized moment line plots...")
    fig3 = figure(3)
    clf()
    plot_3Dsym_standardized(xm, S, S5, nmin, nmax, cc, Np)
    tight_layout()
    if save_figures
        savefig(joinpath(output_dir, "fig03_standardized_moments.png"), dpi=150)
    end
    fig3.canvas.draw()  # Force drawing to ensure figure persists
    
    println("Creating Figure 4: Contour plots...")
    fig4 = figure(4)
    clf()
    contour_plots_3D(xm, ym, M_plot, C, S, Np)
    tight_layout()
    if save_figures
        savefig(joinpath(output_dir, "fig04_contours.png"), dpi=150)
    end
    fig4.canvas.draw()  # Force drawing to ensure figure persists
    
    println("Creating Figure 5: C-moment contour plots...")
    fig5 = figure(5)
    clf()
    Cmoment_plots_3D(xm, ym, S, Np)
    tight_layout()
    if save_figures
        savefig(joinpath(output_dir, "fig05_cmoments.png"), dpi=150)
    end
    fig5.canvas.draw()  # Force drawing to ensure figure persists
    
    println("Creating Figure 6: S-moment contour plots...")
    fig6 = figure(6)
    clf()
    Smoment_plots_3D(xm, ym, S, Np)
    tight_layout()
    if save_figures
        savefig(joinpath(output_dir, "fig06_smoments.png"), dpi=150)
    end
    fig6.canvas.draw()  # Force drawing to ensure figure persists
    
    println("Creating Figure 7: Hyperbolicity plots...")
    fig7 = figure(7)
    clf()
    hyperbolic_plots_3D(xm, ym, eig_data, Np)
    tight_layout()
    if save_figures
        savefig(joinpath(output_dir, "fig07_hyperbolicity.png"), dpi=150)
    end
    fig7.canvas.draw()  # Force drawing to ensure figure persists
    
    println("="^70)
    println("Visualization complete! Created 7 figures.")
    if save_figures
        println("Figures saved to: $output_dir")
        println("  - fig01_moments.png")
        println("  - fig02_central_moments.png")
        println("  - fig03_standardized_moments.png")
        println("  - fig04_contours.png")
        println("  - fig05_cmoments.png")
        println("  - fig06_smoments.png")
        println("  - fig07_hyperbolicity.png")
    else
        println("Figures displayed in windows.")
        println("Press Enter to close all windows and exit...")
        # Show the figures
        PyPlot.show()
        # Wait for user input to keep windows open
        readline()
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
    format_colorbar()
    
    # Panel 2: U velocity
    subplot(3, 4, 2)
    Z = M[:, :, 2] ./ M[:, :, 1]
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("U velocity")
    format_colorbar()
    
    # Panel 3: V velocity
    subplot(3, 4, 3)
    Z = M[:, :, 6] ./ M[:, :, 1]
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("V velocity")
    format_colorbar()
    
    # Panel 4: S_110
    subplot(3, 4, 4)
    Z = S[:, :, 7]
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("S_{110}")
    format_colorbar()
    
    # Panel 5: C_200
    subplot(3, 4, 5)
    Z = C[:, :, 3]
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("C_{200}")
    format_colorbar()
    
    # Panel 6: C_020
    subplot(3, 4, 6)
    Z = C[:, :, 10]
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("C_{020}")
    format_colorbar()
    
    # Panel 7: C_002
    subplot(3, 4, 7)
    Z = C[:, :, 20]
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("C_{002}")
    format_colorbar()
    
    # Panel 8: |Delta_1|
    subplot(3, 4, 8)
    Z = 1 .+ 2 .* S[:, :, 7] .* S[:, :, 17] .* S[:, :, 26] .- 
        S[:, :, 7].^2 .- S[:, :, 17].^2 .- S[:, :, 26].^2
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("|Δ_1|")
    format_colorbar()
    
    # Panel 9: H_200
    subplot(3, 4, 9)
    Z = max.(0, S[:, :, 5] .- S[:, :, 4].^2 .- 1)
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("H_{200}")
    format_colorbar()
    
    # Panel 10: H_020
    subplot(3, 4, 10)
    Z = max.(0, S[:, :, 15] .- S[:, :, 13].^2 .- 1)
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("H_{020}")
    format_colorbar()
    
    # Panel 11: H_002
    subplot(3, 4, 11)
    Z = max.(0, S[:, :, 25] .- S[:, :, 23].^2 .- 1)
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("H_{002}")
    format_colorbar()
    
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
    format_colorbar()
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
        format_colorbar()
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
        format_colorbar()
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
        format_colorbar()
    end
    
    # Plot eigenvalue differences
    subplot(3, 3, 7)
    Z = lam6x[:, :, 2] .- lam6x[:, :, 1]
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("λ_2 - λ_1")
    format_colorbar()
    
    subplot(3, 3, 8)
    Z = lam6x[:, :, 4] .- lam6x[:, :, 3]
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("λ_4 - λ_3")
    format_colorbar()
    
    subplot(3, 3, 9)
    Z = lam6x[:, :, 6] .- lam6x[:, :, 5]
    contourf(X, Y, Z', 50, cmap="sky")
    axis("square")
    title("λ_6 - λ_5")
    format_colorbar()
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
    Y = [M[i, i, 1] for i = 1:Np]
    plot(xm[nmin:nmax], Y[nmin:nmax], color=cc, linewidth=2)
    xlabel("x")
    ylabel("M_{000}")
    grid(true)
    ax = gca()
    ax.tick_params(labelsize=18)
    
    # Panel 2: u_1 (mean velocity)
    subplot(nl, nc, 2)
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
    Y = [S[i, i, 4] for i = 1:Np]
    plot(xm[nmin:nmax], Y[nmin:nmax], color=cc, linewidth=2)
    xlabel("x")
    ylabel("S_{300}")
    grid(true)
    ax = gca()
    ax.tick_params(labelsize=18)
    
    # Panel 2: S_400
    subplot(nl, nc, 2)
    Y = [S[i, i, 5] for i = 1:Np]
    plot(xm[nmin:nmax], Y[nmin:nmax], color=cc, linewidth=2)
    xlabel("x")
    ylabel("S_{400}")
    grid(true)
    ax = gca()
    ax.tick_params(labelsize=18)
    
    # Panel 3: H_400
    subplot(nl, nc, 3)
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
    Y = [S5[i, i, 2] for i = 1:Np]
    plot(xm[nmin:nmax], Y[nmin:nmax], color=cc, linewidth=2)
    xlabel("x")
    ylabel("S_{410}")
    grid(true)
    ax = gca()
    ax.tick_params(labelsize=18)
    
    # Panel 10: S_320
    subplot(nl, nc, 10)
    Y = [S5[i, i, 3] for i = 1:Np]
    plot(xm[nmin:nmax], Y[nmin:nmax], color=cc, linewidth=2)
    xlabel("x")
    ylabel("S_{320}")
    grid(true)
    ax = gca()
    ax.tick_params(labelsize=18)
    
    # Panel 11: |Delta_1|
    subplot(nl, nc, 11)
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
Note: Using 'viridis' as a replacement since 'sky' is not available in all matplotlib versions.
"""
function set_sky_colormap()
    # Create a custom 'sky' colormap that matches MATLAB's sky colormap
    # Sky colormap goes from dark blue -> blue -> light blue -> white
    
    # Import matplotlib.colors for creating custom colormap
    matplotlib = pyimport("matplotlib")
    colors_module = pyimport("matplotlib.colors")
    
    # Check if 'sky' colormap is already registered
    try
        # Try to get the existing colormap
        existing_cmap = matplotlib.colormaps.get_cmap("sky")
        return existing_cmap
    catch
        # If it doesn't exist, create and register it
        # Define the sky colormap colors (RGB values)
        sky_colors = [
            (0.0, 0.0, 0.5),    # Dark blue
            (0.0, 0.3, 0.8),    # Blue  
            (0.3, 0.7, 1.0),    # Light blue
            (0.7, 0.9, 1.0),    # Very light blue
            (1.0, 1.0, 1.0)     # White
        ]
        
        # Create the custom colormap
        sky_cmap = colors_module.LinearSegmentedColormap.from_list("sky", sky_colors, N=256)
        
        # Register the colormap so it can be used by name
        matplotlib.colormaps.register(sky_cmap, name="sky")
        
        return sky_cmap
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

# ============================================================================
# 3D Visualization Functions
# ============================================================================

"""
    plot_multiple_z_slices(M, xm, ym, zm, quantity_name; nslices=6, fignum=nothing)

Create multi-panel plot showing 2D contours at different z-levels.

# Arguments
- `M`: 4D moment array (Np x Np x Nz x Nmom)
- `xm`: X-coordinates of cell centers
- `ym`: Y-coordinates of cell centers  
- `zm`: Z-coordinates of cell centers
- `quantity_name`: Name of quantity to plot ("density", "U", "V", "W", "C200", etc.)
- `nslices`: Number of z-slices to show (default: 6)
- `fignum`: Figure number (optional)

# Returns
- Figure object
"""
function plot_multiple_z_slices(M, xm, ym, zm, quantity_name; nslices=6, fignum=nothing)
    Np = size(M, 1)
    Nz = size(M, 3)
    Nmom = size(M, 4)
    
    # Select z-indices to plot (evenly spaced)
    z_indices = round.(Int, range(1, Nz, length=nslices))
    
    # Create figure
    if fignum !== nothing
        fig = figure(fignum)
    else
        fig = figure()
    end
    clf()
    
    # Determine layout (2 rows for 4-6 slices, 3 rows for more)
    nrows = nslices <= 6 ? 2 : 3
    ncols = ceil(Int, nslices / nrows)
    
    X, Y = meshgrid(xm, ym)
    set_sky_colormap()
    
    for (panel_idx, k) in enumerate(z_indices)
        subplot(nrows, ncols, panel_idx)
        
        # Extract quantity based on name
        if quantity_name == "density"
            Z = M[:, :, k, 1]
            title_str = "Density"
        elseif quantity_name == "U"
            Z = M[:, :, k, 2] ./ M[:, :, k, 1]
            title_str = "U velocity"
        elseif quantity_name == "V"
            Z = M[:, :, k, 6] ./ M[:, :, k, 1]
            title_str = "V velocity"
        elseif quantity_name == "W"
            Z = M[:, :, k, 16] ./ M[:, :, k, 1]  # M001 is at index 16
            title_str = "W velocity"
        elseif quantity_name == "C200"
            # Need to compute C moments
            C = zeros(Np, Np, Nmom)
            for i = 1:Np
                for j = 1:Np
                    MOM = M[i, j, k, :]
                    CC, _ = M2CS4_35(MOM)
                    C[i, j, :] = CC
                end
            end
            Z = C[:, :, 3]
            title_str = "C_{200}"
        elseif quantity_name == "C020"
            C = zeros(Np, Np, Nmom)
            for i = 1:Np
                for j = 1:Np
                    MOM = M[i, j, k, :]
                    CC, _ = M2CS4_35(MOM)
                    C[i, j, :] = CC
                end
            end
            Z = C[:, :, 10]
            title_str = "C_{020}"
        elseif quantity_name == "C002"
            C = zeros(Np, Np, Nmom)
            for i = 1:Np
                for j = 1:Np
                    MOM = M[i, j, k, :]
                    CC, _ = M2CS4_35(MOM)
                    C[i, j, :] = CC
                end
            end
            Z = C[:, :, 20]
            title_str = "C_{002}"
        else
            error("Unknown quantity: $quantity_name")
        end
        
        contourf(X, Y, Z', 50, cmap="sky")
        axis("square")
        title("$(title_str) (z=$(round(zm[k], digits=3)))")
        format_colorbar(shrink=0.9, aspect=15)
        xlabel("x")
        ylabel("y")
    end
    
    tight_layout()
    return fig
end

"""
    plot_3d_isosurface(M, xm, ym, zm, quantity_name, isovalue; fignum=nothing, alpha=0.7)

Create 3D isosurface plot showing surfaces of constant value.

Uses PyPlot's mplot3d toolkit to create true 3D visualization.

# Arguments
- `M`: 4D moment array (Np x Np x Nz x Nmom)
- `xm`: X-coordinates of cell centers
- `ym`: Y-coordinates of cell centers
- `zm`: Z-coordinates of cell centers
- `quantity_name`: Name of quantity ("density", "U", "V", "W")
- `isovalue`: Value at which to draw isosurface
- `fignum`: Figure number (optional)
- `alpha`: Transparency (default: 0.7)

# Returns
- Figure object
"""
function plot_3d_isosurface(M, xm, ym, zm, quantity_name, isovalue; fignum=nothing, alpha=0.7)
    Np = size(M, 1)
    Nz = size(M, 3)
    
    # Import 3D plotting toolkit
    mplot3d = pyimport("mpl_toolkits.mplot3d")
    
    # Create figure
    if fignum !== nothing
        fig = figure(fignum)
    else
        fig = figure(figsize=(10, 8))
    end
    clf()
    
    ax = fig.add_subplot(111, projection="3d")
    
    # Extract quantity
    if quantity_name == "density"
        Q = M[:, :, :, 1]
        title_str = "Density Isosurface"
    elseif quantity_name == "U"
        Q = M[:, :, :, 2] ./ M[:, :, :, 1]
        title_str = "U Velocity Isosurface"
    elseif quantity_name == "V"
        Q = M[:, :, :, 6] ./ M[:, :, :, 1]
        title_str = "V Velocity Isosurface"
    elseif quantity_name == "W"
        Q = M[:, :, :, 16] ./ M[:, :, :, 1]  # M001 is at index 16
        title_str = "W Velocity Isosurface"
    else
        error("Unknown quantity: $quantity_name")
    end
    
    # Use scatter plot for 3D visualization (more reliable than contour3D)
    # Show points where Q ≈ isovalue
    tolerance = max(0.1 * abs(isovalue), 0.01)
    mask = abs.(Q .- isovalue) .< tolerance
    indices = findall(mask)
    
    if !isempty(indices)
        # Extract coordinates
        x_pts = Float64[]
        y_pts = Float64[]
        z_pts = Float64[]
        q_vals = Float64[]
        
        for idx in indices
            i, j, k = idx[1], idx[2], idx[3]
            push!(x_pts, xm[i])
            push!(y_pts, ym[j])
            push!(z_pts, zm[k])
            push!(q_vals, Q[idx])
        end
        
        # Plot as scatter
        scatter = ax.scatter3D(x_pts, y_pts, z_pts, c=q_vals, 
                               alpha=alpha, s=30, cmap="viridis",
                               vmin=isovalue-tolerance, vmax=isovalue+tolerance)
        fig.colorbar(scatter, ax=ax, shrink=0.5, aspect=10)
    else
        @warn "No points found near isovalue=$isovalue for $quantity_name"
    end
    
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title("$(title_str) = $(isovalue)")
    
    return fig
end

"""
    plot_centerline_profiles(M, xm, ym, zm, Np, Nz; fignum=nothing)

Plot line profiles along centerlines showing z-variation.

Shows how density and velocities vary along z at the center of the domain.

# Arguments
- `M`: 4D moment array (Np x Np x Nz x Nmom)
- `xm`: X-coordinates
- `ym`: Y-coordinates
- `zm`: Z-coordinates
- `Np`: Grid size in x-y
- `Nz`: Grid size in z
- `fignum`: Figure number (optional)

# Returns
- Figure object
"""
function plot_centerline_profiles(M, xm, ym, zm, Np, Nz; fignum=nothing)
    # Find center indices
    i_center = div(Np, 2) + 1
    j_center = div(Np, 2) + 1
    
    # Create figure
    if fignum !== nothing
        fig = figure(fignum)
    else
        fig = figure(figsize=(12, 8))
    end
    clf()
    
    # Extract profiles at center (x=0, y=0)
    rho_profile = M[i_center, j_center, :, 1]
    U_profile = M[i_center, j_center, :, 2] ./ M[i_center, j_center, :, 1]
    V_profile = M[i_center, j_center, :, 6] ./ M[i_center, j_center, :, 1]
    W_profile = M[i_center, j_center, :, 16] ./ M[i_center, j_center, :, 1]  # M001 is at index 16
    
    # Plot density
    subplot(2, 2, 1)
    plot(zm, rho_profile, "b-", linewidth=2)
    xlabel("z")
    ylabel("Density")
    title("Density along centerline (x=0, y=0)")
    grid(true)
    
    # Plot U velocity
    subplot(2, 2, 2)
    plot(zm, U_profile, "r-", linewidth=2)
    xlabel("z")
    ylabel("U velocity")
    title("U velocity along centerline")
    grid(true)
    
    # Plot V velocity
    subplot(2, 2, 3)
    plot(zm, V_profile, "g-", linewidth=2)
    xlabel("z")
    ylabel("V velocity")
    title("V velocity along centerline")
    grid(true)
    
    # Plot W velocity
    subplot(2, 2, 4)
    plot(zm, W_profile, "m-", linewidth=2)
    xlabel("z")
    ylabel("W velocity")
    title("W velocity along centerline")
    grid(true)
    
    tight_layout()
    return fig
end

"""
    plot_3d_diagnostics(M, xm, ym, zm, Np, Nz, Nmom; save_figures=false, output_dir=".")

Comprehensive 3D diagnostic visualization suite.

Creates multiple figures showing:
- Multi-slice contours for density and velocities
- C-moment slices
- Realizability metrics (Delta1, Delta2, H200, H020)
- Centerline profiles

# Arguments
- `M`: 4D moment array (Np x Np x Nz x Nmom)
- `xm`: X-coordinates
- `ym`: Y-coordinates
- `zm`: Z-coordinates
- `Np`: Grid size
- `Nz`: Grid size in z
- `Nmom`: Number of moments
- `save_figures`: Save to disk (default: false)
- `output_dir`: Output directory (default: ".")
"""
function plot_3d_diagnostics(M, xm, ym, zm, Np, Nz, Nmom; save_figures=false, output_dir=".")
    println("\n" * "="^70)
    println("Creating 3D Diagnostic Plots...")
    println("="^70)
    
    # Figure 101: Density slices
    println("Figure 101: Density z-slices...")
    fig101 = plot_multiple_z_slices(M, xm, ym, zm, "density", nslices=6, fignum=101)
    fig101.canvas.draw()
    if save_figures
        savefig(joinpath(output_dir, "fig101_density_slices.png"), dpi=150)
    end
    
    # Figure 102: U velocity slices
    println("Figure 102: U velocity z-slices...")
    fig102 = plot_multiple_z_slices(M, xm, ym, zm, "U", nslices=6, fignum=102)
    fig102.canvas.draw()
    if save_figures
        savefig(joinpath(output_dir, "fig102_U_slices.png"), dpi=150)
    end
    
    # Figure 103: V velocity slices
    println("Figure 103: V velocity z-slices...")
    fig103 = plot_multiple_z_slices(M, xm, ym, zm, "V", nslices=6, fignum=103)
    fig103.canvas.draw()
    if save_figures
        savefig(joinpath(output_dir, "fig103_V_slices.png"), dpi=150)
    end
    
    # Figure 104: W velocity slices
    println("Figure 104: W velocity z-slices...")
    fig104 = plot_multiple_z_slices(M, xm, ym, zm, "W", nslices=6, fignum=104)
    fig104.canvas.draw()
    if save_figures
        savefig(joinpath(output_dir, "fig104_W_slices.png"), dpi=150)
    end
    
    # Figure 105: C200 slices
    println("Figure 105: C200 z-slices...")
    fig105 = plot_multiple_z_slices(M, xm, ym, zm, "C200", nslices=6, fignum=105)
    fig105.canvas.draw()
    if save_figures
        savefig(joinpath(output_dir, "fig105_C200_slices.png"), dpi=150)
    end
    
    # Figure 106: C020 slices
    println("Figure 106: C020 z-slices...")
    fig106 = plot_multiple_z_slices(M, xm, ym, zm, "C020", nslices=6, fignum=106)
    fig106.canvas.draw()
    if save_figures
        savefig(joinpath(output_dir, "fig106_C020_slices.png"), dpi=150)
    end
    
    # Figure 107: C002 slices
    println("Figure 107: C002 z-slices...")
    fig107 = plot_multiple_z_slices(M, xm, ym, zm, "C002", nslices=6, fignum=107)
    fig107.canvas.draw()
    if save_figures
        savefig(joinpath(output_dir, "fig107_C002_slices.png"), dpi=150)
    end
    
    # Figure 108: Centerline profiles
    println("Figure 108: Centerline profiles...")
    fig108 = plot_centerline_profiles(M, xm, ym, zm, Np, Nz, fignum=108)
    fig108.canvas.draw()
    if save_figures
        savefig(joinpath(output_dir, "fig108_centerline_profiles.png"), dpi=150)
    end
    
    # Figure 109: Realizability metrics at multiple z-slices
    println("Figure 109: Delta1 and H200 z-slices...")
    fig109 = figure(109)
    clf()
    
    # Select 4 z-slices
    z_indices = round.(Int, range(1, Nz, length=4))
    X, Y = meshgrid(xm, ym)
    set_sky_colormap()
    
    for (panel_idx, k) in enumerate(z_indices)
        # Compute S moments for this slice
        S = zeros(Np, Np, Nmom)
        for i = 1:Np
            for j = 1:Np
                MOM = M[i, j, k, :]
                _, SS = M2CS4_35(MOM)
                S[i, j, :] = SS
            end
        end
        
        # Top row: Delta1
        subplot(2, 4, panel_idx)
        Delta1 = 1 .+ 2 .* S[:, :, 7] .* S[:, :, 17] .* S[:, :, 26] .- 
                 S[:, :, 7].^2 .- S[:, :, 17].^2 .- S[:, :, 26].^2
        contourf(X, Y, Delta1', 50, cmap="sky")
        axis("square")
        title("Δ₁ (z=$(round(zm[k], digits=3)))")
        format_colorbar(shrink=0.9, aspect=15)
        
        # Bottom row: H200
        subplot(2, 4, panel_idx + 4)
        H200 = max.(0, S[:, :, 5] .- S[:, :, 4].^2 .- 1)
        contourf(X, Y, H200', 50, cmap="sky")
        axis("square")
        title("H₂₀₀ (z=$(round(zm[k], digits=3)))")
        format_colorbar(shrink=0.9, aspect=15)
    end
    
    tight_layout()
    fig109.canvas.draw()
    if save_figures
        savefig(joinpath(output_dir, "fig109_realizability_slices.png"), dpi=150)
    end
    
    # Figure 110: 3D Isosurface of density
    println("Figure 110: 3D density isosurface...")
    try
        # Choose isovalue as midpoint between background and jet density
        rho_min = minimum(M[:, :, :, 1])
        rho_max = maximum(M[:, :, :, 1])
        isovalue = 0.3 * (rho_max - rho_min) + rho_min
        
        fig110 = plot_3d_isosurface(M, xm, ym, zm, "density", isovalue, fignum=110, alpha=0.6)
        fig110.canvas.draw()
        if save_figures
            savefig(joinpath(output_dir, "fig110_density_isosurface.png"), dpi=150)
        end
    catch e
        @warn "Could not create 3D isosurface plot" exception=e
    end
    
    println("="^70)
    println("3D diagnostic plots complete!")
    println("="^70)
end

