# This is still in cosntruction so expect bugs.

import math
import os

# Import your existing tools
from rna_loader import load_rna_structure
from rna_distance import get_all_distances

def get_pair_name(res1, res2):
    """
    Standardizes pair names so order doesn't matter.
    Example: 'A' and 'U' -> 'AU', 'U' and 'A' -> 'AU'.
    This reduces the combinations to the required 10 types.
    """
    # Sort the names alphabetically to ensure consistency
    sorted_pair = sorted([res1, res2])
    return f"{sorted_pair[0]}{sorted_pair[1]}"

def train_objective_function(pdb_files, atom_type="C3'", bin_size=1.0, max_dist=20.0):
    """
    Main training logic.
    1. Reads PDB files.
    2. Builds histograms for specific pairs and the reference (XX).
    3. Computes the scoring function.
    """
    
    # 1. Initialize Histograms
    # We need 10 specific pairs + 1 Reference "XX"
    # Each histogram is a list of counts, initialized to 0.
    num_bins = int(max_dist / bin_size) + 1
    
    # Dictionary to store counts: {'AU': [0,0...], 'GC': [0,0...] ...}
    pair_counts = {} 
    ref_counts = [0] * num_bins  # This is the 'XX' histogram
    
    valid_bases = ['A', 'U', 'C', 'G']
    
    # Pre-fill dictionary keys for the 10 canonical pairs
    for i in range(len(valid_bases)):
        for j in range(i, len(valid_bases)): # Start at i to include AA, UU etc.
            pair_name = get_pair_name(valid_bases[i], valid_bases[j])
            pair_counts[pair_name] = [0] * num_bins

    print(f"Initialized training for {len(pdb_files)} files...")

    # 2. Process Files and Count
    for filepath in pdb_files:
        if not os.path.exists(filepath):
            continue
            
        structure = load_rna_structure(filepath)
        if structure is None:
            continue
            
        # Get distances using your existing code
        # Note: Training only uses Intrachain distances 
        interactions = get_all_distances(structure, atom_name=atom_type, max_distance=max_dist)
        
        for item in interactions:
            # Training Requirement: Only Intrachain 
            if item['Type'] != 'Intrachain':
                continue
                
            r1 = item['Res1']
            r2 = item['Res2']
            dist = item['Distance']
            
            # Skip non-standard bases (e.g., modified bases not in A,U,C,G)
            if r1 not in valid_bases or r2 not in valid_bases:
                continue

            # Determine Bin Index
            # e.g., distance 5.9 -> index 5
            bin_idx = int(dist // bin_size)
            
            if bin_idx >= num_bins:
                continue

            # Increment Pair Count (Observed)
            pair_name = get_pair_name(r1, r2)
            pair_counts[pair_name][bin_idx] += 1
            
            # Increment Reference Count (XX) [cite: 102]
            ref_counts[bin_idx] += 1

    # 3. Calculate Scores (Pseudo-Energy)
    # Score = -log( P(obs) / P(ref) )
    print("Computing scores...")
    
    # We return a dictionary containing the final scores for each pair
    # Format: {'AU': [score_bin_0, score_bin_1...], ...}
    final_scores = {}
    
    # Calculate Total Counts (N_ij and N_xx) for normalization
    total_ref = sum(ref_counts)
    
    for pair, counts in pair_counts.items():
        scores = []
        total_pair = sum(counts)
        
        for r in range(num_bins):
            # Avoid division by zero
            if total_pair == 0 or total_ref == 0:
                scores.append(10.0) # Max penalty
                continue
                
            # Frequency observed: count / total_pair
            freq_obs = counts[r] / total_pair
            
            # Frequency reference: count_ref / total_ref
            freq_ref = ref_counts[r] / total_ref
            
            # Calculate Score
            if freq_obs > 0 and freq_ref > 0:
                ratio = freq_obs / freq_ref
                val = -math.log(ratio)
                scores.append(val)
            else:
                # If we never observe this distance, assign high energy (penalty)
                # Requirement: Max scoring value arbitrarily set to 10 [cite: 112]
                scores.append(10.0) 
                
        final_scores[pair] = scores
        
    return final_scores

def save_scores(scores_dict, output_dir="potentials"):
    """
    Writes the 10 files required by instructions.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    for pair, values in scores_dict.items():
        filename = os.path.join(output_dir, f"{pair}.txt")
        with open(filename, "w") as f:
            for v in values:
                f.write(f"{v:.4f}\n")
    print(f"Saved 10 scoring files to folder '{output_dir}'")

# --- usage example ---
if __name__ == "__main__":
    # List your training PDBs here
    # For testing, we just use the one you have
    training_files = ["1EHZ.pdb", "1U9S.pdb"] 
    
    scores = train_objective_function(training_files)
    save_scores(scores)