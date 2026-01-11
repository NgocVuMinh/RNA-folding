# rna_distance.py
# All-atom version - computes distances for ALL atoms in RNA residues
from Bio.PDB import Entity

# Standard RNA nucleotides
VALID_RNA_BASES = {'A', 'U', 'C', 'G'}

def is_strict_pure_rna(chain):
    """
    Checks if a chain contains ONLY standard RNA bases (A, U, C, G).
    Returns False if any protein, DNA, or modified base is found.
    """
    for residue in chain:
        # Ignore water molecules and ions (heteroatoms)
        if residue.id[0] != ' ':
            continue
            
        res_name = residue.get_resname().strip()
        
        # If we find even one non-standard base, reject the whole chain
        if res_name not in VALID_RNA_BASES:
            return False
            
    return True

def get_all_distances(model, atom_name=None, max_distance=20.0, seq_sep=3):
    """
    Calculates pairwise distances for ALL ATOMS in clean RNA chains.
    
    Args:
        model: Bio.PDB Model object.
        atom_name (str): Ignored in all-atom version (for compatibility). 
                         All atoms will be used.
        max_distance (float): Cutoff for interaction (default 20A).
        seq_sep (int): Sequence separation cutoff. 
                            If 3, discards i to i+3, keeps i to i+4.
    """
    interactions = []
    
    # 1. Extract all atoms from pure RNA chains
    valid_atoms = []
    
    for chain in model:
        # Filter: Skip "dirty" chains (containing proteins/modified bases)
        if not is_strict_pure_rna(chain):
            continue

        for residue in chain:
            if residue.id[0] != ' ':
                continue
            
            # Extract ALL atoms in this residue
            for atom in residue:
                valid_atoms.append({
                    "chain": chain.id,
                    "res_num": residue.id[1],  # Sequence number
                    "res_name": residue.get_resname().strip(),
                    "atom_name": atom.name,
                    "atom": atom
                })

    # 2. Compute distances
    count = len(valid_atoms)
    for i in range(count):
        for j in range(i + seq_sep, count):
            atom_A = valid_atoms[i]
            atom_B = valid_atoms[j]
            
            # Only consider intrachain interactions
            if atom_A['chain'] != atom_B['chain']:
                continue
            
            interaction_type = "Intrachain"

            # --- CALCULATION: Euclidean distance ---
            dist = atom_A['atom'] - atom_B['atom']
            
            if dist <= max_distance:
                interactions.append({
                    "Res1": atom_A['res_name'],
                    "Res2": atom_B['res_name'],
                    "AtomA": atom_A['atom_name'],
                    "AtomB": atom_B['atom_name'],
                    "Type": interaction_type,
                    "Distance": dist
                })
                
    return interactions
