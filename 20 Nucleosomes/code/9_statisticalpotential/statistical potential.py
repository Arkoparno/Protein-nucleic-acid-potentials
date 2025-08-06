import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

#  CONFIGURATION 
INPUT_CSV = "interaction_orientations.csv"
OUT_DIR = "statistical_potential_3D"
os.makedirs(OUT_DIR, exist_ok=True)

# Constants
kB_T = 0.593  # kcal/mol at 300K

#  LOAD DATA 
df = pd.read_csv(INPUT_CSV)
r     = df["distance"].values
theta = np.radians(df["theta_deg"].values)
phi   = np.radians(df["phi_deg"].values)

#  DEFINE BINS 
r_bins     = np.linspace(0, 100, 51)        # 0–100 Å, 2Å bins
theta_bins = np.linspace(0, np.pi, 37)      # 0–π rad, 5° bins
phi_bins   = np.linspace(0, 2*np.pi, 73)    # 0–2π rad, 5° bins

r_centers     = (r_bins[:-1] + r_bins[1:]) / 2
theta_centers = (theta_bins[:-1] + theta_bins[1:]) / 2
phi_centers   = (phi_bins[:-1] + phi_bins[1:]) / 2

#  BUILD OBSERVED DISTRIBUTION 
H_obs, _ = np.histogramdd(
    sample=np.vstack((r, theta, phi)).T,
    bins=(r_bins, theta_bins, phi_bins)
)

p_obs = H_obs / np.sum(H_obs)  # Normalize to total probability

#  BUILD REFERENCE DISTRIBUTION 
R, T, P = np.meshgrid(r_centers, theta_centers, phi_centers, indexing='ij')
p_ref = R**2 * np.sin(T)
p_ref /= np.sum(p_ref)  # Normalize to total 1

#  COMPUTE POTENTIAL 
with np.errstate(divide='ignore', invalid='ignore'):
    ratio = np.where((p_obs > 0) & (p_ref > 0), p_obs / p_ref, 1e-8)
    U = -kB_T * np.log(ratio)

#  AVERAGE POTENTIALS 
U_r_avg     = np.nanmean(U, axis=(1, 2))  # average over θ and φ
U_theta_avg = np.nanmean(U, axis=(0, 2))  # average over r and φ
U_phi_avg   = np.nanmean(U, axis=(0, 1))  # average over r and θ

#  SAVE & PLOT RESULTS 
# 1. U(r)
plt.figure(figsize=(8, 5))
plt.plot(r_centers, U_r_avg, color="crimson")
plt.xlabel("Distance r (Å)")
plt.ylabel("U(r) [kcal/mol]")
plt.title("Potential of Mean Force: U(r)")
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/U_r.png", dpi=300)
plt.close()

# 2. U(θ)
plt.figure(figsize=(8, 5))
plt.plot(np.degrees(theta_centers), U_theta_avg, color="darkblue")
plt.xlabel("Polar angle θ (degrees)")
plt.ylabel("U(θ) [kcal/mol]")
plt.title("Potential of Mean Force: U(θ)")
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/U_theta.png", dpi=300)
plt.close()

# 3. U(φ)
plt.figure(figsize=(8, 5))
plt.plot(np.degrees(phi_centers), U_phi_avg, color="seagreen")
plt.xlabel("Azimuthal angle φ (degrees)")
plt.ylabel("U(φ) [kcal/mol]")
plt.title("Potential of Mean Force: U(φ)")
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/U_phi.png", dpi=300)
plt.close()

# Save raw potential values
pd.DataFrame({"r": r_centers, "U(r)": U_r_avg}).to_csv(f"{OUT_DIR}/U_r.csv", index=False)
pd.DataFrame({"theta_deg": np.degrees(theta_centers), "U(theta)": U_theta_avg}).to_csv(f"{OUT_DIR}/U_theta.csv", index=False)
pd.DataFrame({"phi_deg": np.degrees(phi_centers), "U(phi)": U_phi_avg}).to_csv(f"{OUT_DIR}/U_phi.csv", index=False)

print(f"[ok] Statistical potential plots saved to: {OUT_DIR}")
