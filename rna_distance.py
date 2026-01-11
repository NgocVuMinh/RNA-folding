# rna_distance.py
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

def get_all_distances(model, atom_name="C3'", max_distance=20.0, seq_sep=3):
    """
    Calculates pairwise distances for clean RNA chains.
    
    Args:
        model: Bio.PDB Model object.
        atom_name (str): Atom to measure (default C3').
        max_distance (float): Cutoff for interaction (default 20A).
        seq_sep (int): Sequence separation cutoff. 
                            If 3, discards i to i+3, keeps i to i+4.
    """
    interactions = []
    
    # 1. Extract valid atoms from pure RNA chains
    valid_atoms = []
    
    for chain in model:
        # Filter: Skip "dirty" chains (containing proteins/modified bases)
        if not is_strict_pure_rna(chain):
            continue

        for residue in chain:
            if residue.id[0] != ' ':
                continue
            
            if atom_name in residue:
                atom = residue[atom_name]
                valid_atoms.append({
                    "chain": chain.id,
                    "res_num": residue.id[1], # Sequence number
                    "res_name": residue.get_resname().strip(),
                    "atom": atom
                })

    # 2. Compute distances
    count = len(valid_atoms)
    for i in range(count):
        for j in range(i + seq_sep, count): # skip 3 neighbors
            atom_A = valid_atoms[i]
            atom_B = valid_atoms[j]
            
            # --- FILTER: Sequence Separation ---
            # if atom_A['chain'] == atom_B['chain']:
            #     seq_dist = abs(atom_A['res_num'] - atom_B['res_num'])
                
            #     # Instruction: "separated by at least 3 positions"
            #     # Logic: If min_dist=3, we skip 1-2, 1-3, 1-4. We keep 1-5 (dist 4).
            #     if seq_dist <= seq_sep:
            #         continue
            #     interaction_type = "Intrachain"
            # else:
            #     interaction_type = "Interchain"
            if atom_A['chain'] != atom_B['chain']:
                continue
            interaction_type = "Intrachain"

            # --- CALCULATION: Euclidean distance ---
            dist = atom_A['atom'] - atom_B['atom']
            
            if dist <= max_distance:
                interactions.append({
                    "Res1": atom_A['res_name'],
                    "Res2": atom_B['res_name'],
                    "Type": interaction_type,
                    "Distance": dist
                })
                
    return interactions