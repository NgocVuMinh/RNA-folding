# rna_distance.py

from Bio.PDB import Entity

# Define the standard RNA alphabet
VALID_RNA_BASES = {'A', 'U', 'C', 'G'}

def is_strict_pure_rna(chain):
    """
    STRICT FILTER: Returns True ONLY if every single residue in the chain is A, U, C, or G.
    
    If it finds a protein residue (ALA, GLY), DNA (DA, DT), or modified RNA (5MC, PSU, H2U), it returns False.
    """
    for residue in chain:
        # Ignore water/ions (heteroatoms) - we only check the main polymer
        if residue.id[0] != ' ':
            continue
            
        res_name = residue.get_resname().strip()
        
        # If we find ONE bad residue, the whole chain is "contaminated"
        if res_name not in VALID_RNA_BASES:
            # print(f"Chain {chain.id} rejected: contains {res_name}")
            return False
            
    return True

def get_all_distances(model, atom_name="C3'", max_distance=20.0, min_distance=3):
    """
    Calculates distances but strictly discards any chain containing 
    non-standard residues.
    """
    interactions = []
    
    # 1. Filter and extract atoms
    valid_atoms = []
    
    for chain in model:
        # --- STRICT CHAIN FILTER ---
        # If the chain contains ANYTHING except A, U, C, G -> skip
        if not is_strict_pure_rna(chain):
            continue

        # If we passed the check, extract all C3' atoms (default) from this pure chain
        for residue in chain:
            if residue.id[0] != ' ':
                continue
            
            if atom_name in residue:
                atom = residue[atom_name]
                valid_atoms.append({
                    "chain": chain.id,
                    "res_num": residue.id[1],
                    "res_name": residue.get_resname().strip(),
                    "atom": atom
                })

    # 2. Calculate pairwise distances
    count = len(valid_atoms)
    for i in range(count):
        for j in range(i + 1, count):
            atom_A = valid_atoms[i]
            atom_B = valid_atoms[j]
            
            # Sequence separation filter (intrachain only)
            is_same_chain = (atom_A['chain'] == atom_B['chain'])
            
            if is_same_chain:
                seq_dist = abs(atom_A['res_num'] - atom_B['res_num'])
                if seq_dist <= min_distance:
                    continue
                interaction_type = "Intrachain"
            else:
                interaction_type = "Interchain"

            # Euclidean distance
            dist = atom_A['atom'] - atom_B['atom']
            
            if dist <= max_distance:
                interactions.append({
                    "Res1": atom_A['res_name'],
                    "Res2": atom_B['res_name'],
                    "Chain1": atom_A['chain'],
                    "Chain2": atom_B['chain'],
                    "Type": interaction_type,
                    "Distance": dist
                })
                
    return interactions
