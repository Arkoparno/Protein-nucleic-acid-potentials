import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# CONFIG 
COM_DIR        = "pdb_separate_COMs"
OUTPUT_CSV     = "interactions.csv"
PLOT_PATH      = "plots interactions.png"
MIN_VALID_DISTANCE = 0.0  

os.makedirs("plots", exist_ok=True)

# UTILITIES 
def euclidean(a, b):
    return np.linalg.norm(np.array(a) - np.array(b))

#  GATHER INTERACTIONS 
records = []
for fname in os.listdir(COM_DIR):
    if not fname.endswith("_protein_COM.csv"):
        continue
    pdb_id = fname.replace("_protein_COM.csv", "")
    protein_file = os.path.join(COM_DIR, f"{pdb_id}_protein_COM.csv")
    dna_file     = os.path.join(COM_DIR, f"{pdb_id}_dna_COM.csv")

    if not os.path.exists(protein_file) or not os.path.exists(dna_file):
        continue

    prot_df = pd.read_csv(protein_file)
    dna_df  = pd.read_csv(dna_file)

    for _, p in prot_df.iterrows():
        bb = (p.bb_x, p.bb_y, p.bb_z)
        sc = (p.sc_x, p.sc_y, p.sc_z)

        for _, d in dna_df.iterrows():
            P = (d.p_x, d.p_y, d.p_z)
            S = (d.s_x, d.s_y, d.s_z)
            B = (d.b_x, d.b_y, d.b_z)

            pairs = {
                "bb_to_p": euclidean(bb, P),
                "bb_to_s": euclidean(bb, S),
                "bb_to_b": euclidean(bb, B),
                "sc_to_p": euclidean(sc, P),
                "sc_to_s": euclidean(sc, S),
                "sc_to_b": euclidean(sc, B),
            }

            for itype, dist in pairs.items():
                if dist >= MIN_VALID_DISTANCE:  
                    records.append({
                        "pdb_id": pdb_id,
                        "aa_type": p.res_name, "aa_resid": p.res_id, "aa_chain": p.chain_id,
                        "dna_type": d.res_name, "dna_resid": d.res_id, "dna_chain": d.chain_id,
                        "interaction_type": itype,
                        "distance": round(dist, 3)
                    })

# EXPORT CSV 
df = pd.DataFrame(records)
df.to_csv(OUTPUT_CSV, index=False)
print(f"[✓] Saved cleaned interactions to {OUTPUT_CSV} (Total: {len(df)})")

# HISTOGRAM 
plt.figure(figsize=(9, 5))
plt.hist(df["distance"], bins=80, color="teal", edgecolor="black", alpha=0.75)
plt.xlabel("Distance (Å)")
plt.ylabel("Interaction Count")
plt.title("Filtered Protein–DNA Center-of-Mass Distances")
plt.grid(True, linestyle="--", alpha=0.4)
plt.tight_layout()
plt.savefig(PLOT_PATH, dpi=300)
plt.show()
