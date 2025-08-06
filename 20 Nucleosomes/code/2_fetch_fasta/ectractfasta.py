from Bio.PDB import PDBParser, PPBuilder
from Bio import SeqIO
import os

input_folder = "/Users/arko/Desktop/nucleosome/nucleo_dna_pdbs"  # your folder of PDBs
output_fasta = "/Users/arko/Desktop/nucleosome/all_proteins.fasta"  # output FASTA file

parser = PDBParser(QUIET=True)
ppb = PPBuilder()
fasta_records = []

for filename in os.listdir(input_folder):
    if filename.endswith(".pdb") or filename.endswith(".cif"):
        pdb_id = filename[:4]
        file_path = os.path.join(input_folder, filename)
        try:
            structure = parser.get_structure(pdb_id, file_path)
            for model in structure:
                for chain in model:
                    sequence = ""
                    for pp in ppb.build_peptides(chain):
                        sequence += str(pp.get_sequence())
                    if sequence:
                        chain_id = chain.id
                        fasta_id = f"{pdb_id}_{chain_id}"
                        fasta_records.append(f">" + fasta_id + "\n" + sequence)
        except Exception as e:
            print(f"Error parsing {filename}: {e}")

with open(output_fasta, "w") as f:
    f.write("\n".join(fasta_records))

print(f"Extracted {len(fasta_records)} protein chains into {output_fasta}")
