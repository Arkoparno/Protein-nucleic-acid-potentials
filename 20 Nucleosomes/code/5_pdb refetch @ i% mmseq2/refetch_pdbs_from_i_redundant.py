import os
import shutil
import re

#  CONFIGURATION 
fasta_file = 'non_redundant_30pc.fasta'
all_pdb_dir = 'nucleo_dna_pdbs'           # Folder with all pdb files
output_dir = 'filtered_pdbs_30'       # Folder to store filtered pdbs

#  CREATE OUTPUT DIRECTORY IF NOT EXISTS 
os.makedirs(output_dir, exist_ok=True)

#  EXTRACT PDB IDs FROM FASTA 
pdb_ids = set()
with open(fasta_file, 'r') as file:
    for line in file:
        if line.startswith('>'):
            match = re.match(r'>?(\w{4})_', line)
            if match:
                pdb_ids.add(match.group(1).lower())  # store as lowercase to match filenames

print(f"Found {len(pdb_ids)} unique PDB IDs in FASTA.")

#  COPY MATCHING PDB FILES 
copied = 0
for pdb_id in pdb_ids:
    pdb_filename = f"{pdb_id}.pdb"
    src_path = os.path.join(all_pdb_dir, pdb_filename)
    dst_path = os.path.join(output_dir, pdb_filename)

    if os.path.exists(src_path):
        shutil.copy2(src_path, dst_path)
        copied += 1
    else:
        print(f"[Warning] {pdb_filename} not found in {all_pdb_dir}")

print(f"\n Copied {copied} matching PDB files to {output_dir}")
