import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load your CSV file
df = pd.read_csv("interaction_orientations.csv")

# Use φ (phi) column
phi = df["phi_deg"]

# Binning parameters
bin_width = 10  # degrees
bins = np.arange(0, 360 + bin_width, bin_width)
bin_centers = (bins[:-1] + bins[1:]) / 2

# Histogram of φ
hist, _ = np.histogram(phi, bins=bins)

# Normalize to get g(φ): relative frequency
g_phi = hist / np.sum(hist)

# Plotting g(φ)
plt.figure(figsize=(10, 5))
plt.plot(bin_centers, g_phi, color='indigo', lw=2, marker='o', alpha=0.8)
plt.xlabel("Azimuthal Angle φ (degrees)", fontsize=12)
plt.ylabel("g(φ)", fontsize=12)
plt.title("Angular Distribution Function g(φ)", fontsize=14)
plt.xticks(np.arange(0, 361, 30))
plt.grid(True, linestyle="--", alpha=0.4)
plt.tight_layout()
plt.savefig("g_phi_plot.png", dpi=300)
plt.show()
