import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load CSV
df = pd.read_csv("interactions.csv")

# Parameters
r_max = 150      # Max distance in Å
bin_width = 1.0  # Bin size in Å
bins = np.arange(0, r_max + bin_width, bin_width)
bin_centers = (bins[:-1] + bins[1:]) / 2

# Step 1: Count number of distances in each bin
hist, _ = np.histogram(df["distance"], bins=bins)

# Step 2: Normalize with shell volume and density
# Estimate system volume from max distance cube (crude approximation)
V_total = (r_max)**3  # in Å³

# Number of total pairs
N_total = len(df)

# Average number density
rho = N_total / V_total

# Shell volumes
volumes = (4/3) * np.pi * (bins[1:]**3 - bins[:-1]**3)

# g(r) calculation
g_r = hist / (rho * volumes)

# Step 3: Plot g(r)
plt.figure(figsize=(8, 5))
plt.plot(bin_centers, g_r, color='darkorange', lw=2)
plt.xlabel("Distance r (Å)", fontsize=12)
plt.ylabel("g(r)", fontsize=12)
plt.title("Radial Distribution Function g(r)", fontsize=14)
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig("gr_rdf_plot.png", dpi=300)
plt.show()
