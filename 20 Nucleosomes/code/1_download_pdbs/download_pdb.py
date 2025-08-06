import os
import requests

# Function to check and download valid PDB ID
def download_pdb(pdb_id, save_dir="pdb_files", file_format="pdb"):
    pdb_id = pdb_id.strip().lower()
    if len(pdb_id) != 4:
        print(f"Skipping invalid ID: {pdb_id}")
        return

    # File extension
    ext = ".pdb" if file_format == "pdb" else ".cif"
    
    # URL based on format
    if file_format == "pdb":
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    elif file_format == "mmCif" or file_format == "cif":
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
    else:
        print(f"Unknown file format: {file_format}")
        return
    
    # Create output directory
    os.makedirs(save_dir, exist_ok=True)
    file_path = os.path.join(save_dir, f"{pdb_id}{ext}")

    # Download
    response = requests.get(url)
    if response.status_code == 200:
        with open(file_path, "wb") as f:
            f.write(response.content)
        print(f"Downloaded {pdb_id.upper()} successfully.")
    else:
        print(f"Failed to download {pdb_id.upper()}: Not found (HTTP {response.status_code})")

# Read PDB IDs
with open("pdb_ids.txt") as file:
    pdb_ids = [line.strip() for line in file if line.strip()]

# Download each
for pdb_id in pdb_ids:
    download_pdb(pdb_id, save_dir="pdb_files", file_format="pdb")  # or file_format="cif"
