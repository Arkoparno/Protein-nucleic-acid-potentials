import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic_2d
from matplotlib import colors
from mpl_toolkits.mplot3d import Axes3D

# Constants
kB_T = 0.593  # kcal/mol at 300K

# Bin settings
r_bins = np.linspace(0, 150, 151)
theta_bins = np.linspace(0, 180, 37)
phi_bins = np.linspace(0, 360, 73)

# Create output base directory
base_dir = "aa_potentials_combinedDNA"
os.makedirs(base_dir, exist_ok=True)

# Load combined CSV
combined_df = pd.read_csv("all_pdb_orientations_combined.csv")

# Utility function to compute U(x) from g(x)
def compute_potential(hist, bins, N, volume_func):
    centers = (bins[:-1] + bins[1:]) / 2
    volume_elements = volume_func(centers)
    g = (hist / N) / (volume_elements / volume_elements.sum())
    U = -kB_T * np.log(np.clip(g, 1e-8, None))
    return centers, g, U

# Heatmap plotting function
def plot_log_heatmap(x, y, u, xlabel, ylabel, title, outname, bins=30, cmap='PuBu_r'):
    u_shifted = u.copy()
    min_u = u.min()
    if min_u <= 0:
        u_shifted += abs(min_u) + 1e-3

    stat, xedges, yedges, _ = binned_statistic_2d(x, y, u_shifted, statistic='mean', bins=bins)
    stat = np.nan_to_num(stat, nan=np.nanmedian(u_shifted))

    X, Y = np.meshgrid(xedges[:-1], yedges[:-1])

    fig, ax = plt.subplots(figsize=(7, 6))
    norm = colors.LogNorm(vmin=stat[stat > 0].min(), vmax=stat.max())
    pcm = ax.pcolor(X, Y, stat.T, cmap=cmap, norm=norm, shading='auto')
    cbar = fig.colorbar(pcm, ax=ax, extend='max')
    cbar.set_label("U_total (shifted, kcal/mol)", fontsize=12)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14)
    plt.tight_layout()
    fig.savefig(outname, dpi=300)
    plt.close(fig)

# 3D plotting function
def plot_3d_potential(r, theta, phi, U, outname):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    sc = ax.scatter(r, theta, phi, c=U, cmap='plasma', s=5)
    ax.set_xlabel('r (Å)')
    ax.set_ylabel('θ (°)')
    ax.set_zlabel('φ (°)')
    plt.colorbar(sc, label='U_total (kcal/mol)')
    ax.set_title(f"3D Potential U(r, θ, φ)\n{os.path.basename(outname).replace('_', ' ').replace('.png','')}", fontsize=11)
    plt.tight_layout()
    fig.savefig(outname, dpi=300)
    plt.close(fig)

# Volume functions
shell_vol = lambda r: 4 * np.pi * r**2 * (r_bins[1] - r_bins[0])
solid_angle_theta = lambda theta: np.sin(np.radians(theta)) * np.radians(theta_bins[1] - theta_bins[0])
solid_angle_phi = lambda phi: np.radians(phi_bins[1] - phi_bins[0])

# Process each amino acid type
for aa in combined_df["aa_type"].unique():
    aa_df = combined_df[combined_df["aa_type"] == aa]
    aa_dir = os.path.join(base_dir, aa)
    os.makedirs(aa_dir, exist_ok=True)

    # DNA combined: DA, DT, DG, DC merged
    dna_df = aa_df[aa_df["dna_type"].isin(["DA", "DT", "DG", "DC"])]

    for inter in ["bb_to_P", "bb_to_S", "bb_to_B", "sc_to_P", "sc_to_S", "sc_to_B"]:
        inter_df = dna_df[dna_df["interaction"] == inter]
        if inter_df.empty:
            continue
        inter_dir = os.path.join(aa_dir, inter)
        os.makedirs(inter_dir, exist_ok=True)

        r_data = inter_df["distance"].values
        theta_data = inter_df["theta_deg"].values
        phi_data = inter_df["phi_deg"].values
        N = len(inter_df)

        # Compute g(r), g(theta), g(phi) and potentials
        hist_r, _ = np.histogram(r_data, bins=r_bins)
        r_centers, g_r, U_r = compute_potential(hist_r, r_bins, N, shell_vol)

        hist_theta, _ = np.histogram(theta_data, bins=theta_bins)
        theta_centers, g_theta, U_theta = compute_potential(hist_theta, theta_bins, N, solid_angle_theta)

        hist_phi, _ = np.histogram(phi_data, bins=phi_bins)
        phi_centers, g_phi, U_phi = compute_potential(hist_phi, phi_bins, N, solid_angle_phi)

        # Build 3D grid and compute total U(r, θ, φ)
        R, THETA, PHI = np.meshgrid(r_centers, theta_centers, phi_centers, indexing="ij")
        U_total = U_r[:, None, None] + U_theta[None, :, None] + U_phi[None, None, :]

        # Flatten for 2D/3D plots
        r_flat = R.ravel()
        theta_flat = THETA.ravel()
        phi_flat = PHI.ravel()
        u_flat = U_total.ravel()

        # Save 2D plots with descriptive titles and consistent names
        plot_log_heatmap(r_flat, theta_flat, u_flat,
                         "r (Å)", "θ (°)",
                         f"U(r, θ) | {aa} to DNA (all) via {inter}",
                         os.path.join(inter_dir, f"{aa}_DNA_all_{inter}_heatmap_r_theta.png"))

        plot_log_heatmap(r_flat, phi_flat, u_flat,
                         "r (Å)", "φ (°)",
                         f"U(r, φ) | {aa} to DNA (all) via {inter}",
                         os.path.join(inter_dir, f"{aa}_DNA_all_{inter}_heatmap_r_phi.png"))

        plot_log_heatmap(theta_flat, phi_flat, u_flat,
                         "θ (°)", "φ (°)",
                         f"U(θ, φ) | {aa} to DNA (all) via {inter}",
                         os.path.join(inter_dir, f"{aa}_DNA_all_{inter}_heatmap_theta_phi.png"))

        # Save 3D plot with title
        plot_3d_potential(r_flat, theta_flat, phi_flat, u_flat,
                          os.path.join(inter_dir, f"{aa}_DNA_all_{inter}_3D_U_total_plot.png"))
