import pandas as pd
import matplotlib.pyplot as plt

# Set font for better compatibility
plt.rcParams["font.family"] = "DejaVu Serif"

# Load the data
count_df = pd.read_csv("protein_dna_counts.csv")
group_df = pd.read_csv("protein_dna_groups.csv")

# Merge on Identity Percentage
merged_df = pd.merge(count_df, group_df, left_on="Identity_Percentage", right_on="Identity(%)")

# Plotting
plt.figure(figsize=(12, 6))

# Plot raw FASTA Count
plt.plot(
    merged_df["Identity_Percentage"],
    merged_df["Fasta_Count"],
    label="FASTA Count",
    color="black",
    linewidth=2
)

# Define colors for each amino acid group
colors = {
    "Aliphatic (%)": "purple",
    "Aromatic (%)": "orange",
    "Polar Uncharged (%)": "green",
    "Acidic (%)": "red",
    "Basic (%)": "blue",
    "Histidine (%)": "brown"
}

# Plot actual percentage values
for group, color in colors.items():
    plt.plot(
        merged_df["Identity_Percentage"],
        merged_df[group],
        label=group,
        color=color,
        linestyle='--'
    )

# Labels and formatting
plt.title("Number of FASTA and % Amino Acids vs % Redundancy (Protein-DNA)")
plt.xlabel("Sequence Identity (%)")
plt.ylabel("Raw Values")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

# Save and show the plot
plt.savefig("protein_dna_combined_graph_RAW.png")
plt.show()
