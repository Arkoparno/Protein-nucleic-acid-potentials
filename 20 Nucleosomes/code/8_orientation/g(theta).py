import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load CSV
df = pd.read_csv("interaction_orientations.csv")

# θ range: 0–180°, bin width: 10°
theta_bins = np.arange(0, 190, 10)
theta_centers = (theta_bins[:-1] + theta_bins[1:]) / 2

# Histogram of θ
hist, _ = np.histogram(df["theta_deg"], bins=theta_bins)


# Normalize by sin(θ) to correct spherical sampling bias
sin_theta = np.sin(np.radians(theta_centers))
g_theta = hist / sin_theta

# Normalize g(θ) so its peak is 1 (optional)
g_theta = g_theta / np.max(g_theta)

# Plot g(θ)
plt.figure(figsize=(8, 5))
plt.plot(theta_centers, g_theta, color="mediumblue", marker="o", lw=2)
plt.xlabel("Polar Angle θ (degrees)", fontsize=12)
plt.ylabel("g(θ) (normalized)", fontsize=12)
plt.title("Angular Distribution Function g(θ)", fontsize=14)
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig("g_theta_plot.png", dpi=300)
plt.show()
