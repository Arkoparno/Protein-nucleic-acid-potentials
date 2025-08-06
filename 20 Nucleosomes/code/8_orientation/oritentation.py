import os
import csv
import math
import numpy as np
import pandas as pd

# ─── CONFIGURATION ──────────────────────────────────────────────
COM_DIR            = "pdb_separate_COMs"
INTERACTIONS_CSV   = "interactions.csv"
OUTPUT_CSV         = "interaction_orientations.csv"
MIN_VALID_DISTANCE = 0.0  # Å — keep all vectors, even close ones

# ─── TRACK ALREADY WRITTEN PDBs ─────────────────────────────────
already_done = set()
if os.path.exists(OUTPUT_CSV):
    with open(OUTPUT_CSV, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            already_done.add(row["pdb_id"])

# ─── CREATE HEADER IF FILE DOESN’T EXIST ────────────────────────
if not os.path.exists(OUTPUT_CSV):
    with open(OUTPUT_CSV, "w", newline="") as outf:
        writer = csv.writer(outf)
        writer.writerow([
            "pdb_id",
            "aa_type", "aa_resid", "aa_chain",
            "dna_type", "dna_resid", "dna_chain",
            "interaction_type", "distance",
            "theta_deg", "phi_deg"
        ])

# ─── ORIENTATION COMPUTATION ───────────────────────────────────
def compute_orientations_for_block(pdb_id, block_rows):
    prot_path = os.path.join(COM_DIR, f"{pdb_id}_protein_COM.csv")
    dna_path  = os.path.join(COM_DIR, f"{pdb_id}_dna_COM.csv")

    if not os.path.exists(prot_path) or not os.path.exists(dna_path):
        print(f"[!] Missing COM files for {pdb_id}, skipping.")
        return []

    prot_df = pd.read_csv(prot_path, dtype={"res_id": str})
    dna_df  = pd.read_csv(dna_path,  dtype={"res_id": str})

    orientations = []

    for row in block_rows:
        aa_resid   = row["aa_resid"]
        aa_chain   = row["aa_chain"]
        dna_resid  = row["dna_resid"]
        dna_chain  = row["dna_chain"]
        itype      = row["interaction_type"]
        dist       = float(row["distance"])

        prot_match = prot_df[(prot_df["res_id"] == aa_resid) & (prot_df["chain_id"] == aa_chain)]
        dna_match  = dna_df[(dna_df["res_id"] == dna_resid) & (dna_df["chain_id"] == dna_chain)]

        if prot_match.empty or dna_match.empty:
            continue

        prot_row = prot_match.iloc[0]
        dna_row  = dna_match.iloc[0]

        try:
            prot_part = np.array([
                prot_row[f"{itype.split('_to_')[0].lower()}_x"],
                prot_row[f"{itype.split('_to_')[0].lower()}_y"],
                prot_row[f"{itype.split('_to_')[0].lower()}_z"]
            ])
            dna_part = np.array([
                dna_row[f"{itype.split('_to_')[1].lower()}_x"],
                dna_row[f"{itype.split('_to_')[1].lower()}_y"],
                dna_row[f"{itype.split('_to_')[1].lower()}_z"]
            ])
        except KeyError as e:
            print(f"[!] {pdb_id} missing expected column for {itype}: {e}")
            continue

        vec = dna_part - prot_part
        r = np.linalg.norm(vec)

         

        if r == 0 or not np.isfinite(r):
            continue  # Skip zero or invalid vectors

        theta = math.degrees(math.acos(np.clip(vec[2] / r, -1.0, 1.0)))
        phi   = math.degrees(math.atan2(vec[1], vec[0])) % 360

        

        orientations.append([
            pdb_id,
            row["aa_type"], aa_resid, aa_chain,
            row["dna_type"], dna_resid, dna_chain,
            itype, dist,
            round(theta, 3), round(phi, 3)
        ])

    return orientations

# ─── MAIN SCRIPT ────────────────────────────────────────────────
with open(INTERACTIONS_CSV, newline="") as inf, \
     open(OUTPUT_CSV, "a", newline="") as outf:

    reader = csv.DictReader(inf)
    writer = csv.writer(outf)

    current_pdb = None
    buffer = []

    for row in reader:
        pdb_id = row["pdb_id"]

        if pdb_id in already_done:
            continue

        if current_pdb is None:
            current_pdb = pdb_id

        if pdb_id != current_pdb:
            results = compute_orientations_for_block(current_pdb, buffer)
            for r in results:
                writer.writerow(r)
            print(f"[✓] {current_pdb}: wrote {len(results)} orientations")
            already_done.add(current_pdb)
            buffer = []
            current_pdb = pdb_id

        buffer.append(row)

    # Final block
    if buffer and current_pdb not in already_done:
        results = compute_orientations_for_block(current_pdb, buffer)
        for r in results:
            writer.writerow(r)
        print(f"[✓] {current_pdb}: wrote {len(results)} orientations")

print("\n✅ All done — orientations saved to:", OUTPUT_CSV)
