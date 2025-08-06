import os
import csv
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

# ========== USER SETTINGS ==========
PDB_FOLDER    = "pdbs"               # <-- point this at your folder of .pdb files
OUTPUT_FOLDER = "pdb_separate_COMs"
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# ========= ATOM CATEGORIES =========
# only standard deoxy‑nucleotides
DNA_RES_NAMES = {"DA", "DT", "DG", "DC"}
DNA_ATOMS = {
    'P': ['P'],
    'S': ["C1'", "C2'", "C3'", "C4'", "C5'", "O4'", "O3'"],
    'B': ["N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"]
}

BACKBONE_ATOMS = {'N', 'CA', 'C', 'O'}
# -----------------------------------

def get_com(atom_list):
    """Center of mass (unweighted) of a list of Atom objects."""
    if not atom_list:
        return (0.0, 0.0, 0.0)
    coords = [atom.get_coord() for atom in atom_list]
    xs, ys, zs = zip(*coords)
    return (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs))

def extract_dna_coms(model):
    rows = []
    for chain in model:
        for residue in chain:
            if residue.resname not in DNA_RES_NAMES:
                continue
            resname  = residue.resname
            resid    = residue.id[1]
            chain_id = chain.id

            groups = {'P': [], 'S': [], 'B': []}
            for atom in residue:
                name = atom.get_name()
                if   name in DNA_ATOMS['P']:
                    groups['P'].append(atom)
                elif name in DNA_ATOMS['S']:
                    groups['S'].append(atom)
                elif name in DNA_ATOMS['B']:
                    groups['B'].append(atom)

            p_com = get_com(groups['P'])
            s_com = get_com(groups['S'])
            b_com = get_com(groups['B'])
            rows.append([
                resname, resid, chain_id,
                p_com[0], p_com[1], p_com[2],
                s_com[0], s_com[1], s_com[2],
                b_com[0], b_com[1], b_com[2],
            ])
    return rows

def extract_protein_coms(model):
    rows = []
    for chain in model:
        for residue in chain:
            # only standard amino acids
            if residue.id[0] != " " or not is_aa(residue):
                continue
            resname  = residue.resname
            resid    = residue.id[1]
            chain_id = chain.id

            bb_atoms = []
            sc_atoms = []
            for atom in residue:
                if atom.get_name() in BACKBONE_ATOMS:
                    bb_atoms.append(atom)
                else:
                    sc_atoms.append(atom)

            bb_com = get_com(bb_atoms)
            sc_com = get_com(sc_atoms)
            rows.append([
                resname, resid, chain_id,
                bb_com[0], bb_com[1], bb_com[2],
                sc_com[0], sc_com[1], sc_com[2],
            ])
    return rows

def write_csv(path, header, data):
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(data)

# ========== MAIN ==========
parser = PDBParser(QUIET=True)

for fn in os.listdir(PDB_FOLDER):
    if not fn.lower().endswith(".pdb"):
        continue

    pdb_id = os.path.splitext(fn)[0]
    in_path = os.path.join(PDB_FOLDER, fn)

    try:
        structure = parser.get_structure(pdb_id, in_path)
        model     = structure[0]

        # DNA/RNA COMs
        dna_rows = extract_dna_coms(model)
        if dna_rows:
            write_csv(
                os.path.join(OUTPUT_FOLDER, f"{pdb_id}_dna_COM.csv"),
                ['res_name','res_id','chain_id',
                 'p_x','p_y','p_z',
                 's_x','s_y','s_z',
                 'b_x','b_y','b_z'],
                dna_rows
            )

        # Protein COMs
        prot_rows = extract_protein_coms(model)
        if prot_rows:
            write_csv(
                os.path.join(OUTPUT_FOLDER, f"{pdb_id}_protein_COM.csv"),
                ['res_name','res_id','chain_id',
                 'bb_x','bb_y','bb_z',
                 'sc_x','sc_y','sc_z'],
                prot_rows
            )

        print(f"✅ {pdb_id} done")

    except Exception as e:
        print(f"❌ {pdb_id} error: {e}")
