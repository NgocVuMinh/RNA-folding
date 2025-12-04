# RNA Training Script
# This script trains an RNA objective function based on provided PDB or CIF files.
# Run this script from the command line with the following options:
# -d <data_folder>: Path to the folder containing structure files.  
# -f <file_format>: The file format to use for training (pdb or cif). Default is pdb.
# Example usage:
# py rna_training.py -d rna_data/pdb -f pdb

import math
import os
import glob
import argparse
import sys

from rna_loader import load_rna_structure
from rna_distance import get_all_distances


def get_pair_name(res1, res2):
    """
    Standardizes pair names so order doesn't matter.
    Example: 'A' and 'U' -> 'AU'.
    """
    sorted_pair = sorted([res1, res2])
    return f"{sorted_pair[0]}{sorted_pair[1]}"


def train_objective_function(structure_files, atom_type="C3'", mode="histogram", bin_size=1.0, min_dist=3.0, max_dist=20.0):
    """
    Main training logic.
    """
    # 1. Initialize histograms
    num_bins = int(max_dist / bin_size) + 1
    pair_counts = {} 
    ref_counts = [0] * num_bins  
    valid_bases = ['A', 'U', 'C', 'G']
    
    for i in range(len(valid_bases)):
        for j in range(i, len(valid_bases)): 
            pair_name = get_pair_name(valid_bases[i], valid_bases[j])
            pair_counts[pair_name] = [0] * num_bins

    print(f"Initialized training for {len(structure_files)} files...")

    # 2. Process files
    for filepath in structure_files:
        if not os.path.exists(filepath):
            continue
            
        structure = load_rna_structure(filepath)
        if structure is None:
            continue
            
        interactions = get_all_distances(structure, atom_name=atom_type, max_distance=max_dist)
        
        for item in interactions:
            if item['Type'] != 'Intrachain':
                continue
                
            r1 = item['Res1']
            r2 = item['Res2']
            dist = item['Distance']
            
            if r1 not in valid_bases or r2 not in valid_bases:
                continue

            bin_idx = int(dist // bin_size)
            if bin_idx >= num_bins:
                continue

            pair_name = get_pair_name(r1, r2)
            pair_counts[pair_name][bin_idx] += 1
            ref_counts[bin_idx] += 1 # [cite: 30]

    # 3. Calculate scores
    print("Computing scores...")
    final_scores = {}
    total_ref = sum(ref_counts)
    
    for pair, counts in pair_counts.items():
        scores = []
        total_pair = sum(counts)
        
        for r in range(num_bins):
            if total_pair == 0 or total_ref == 0:
                scores.append(10.0) 
                continue
                
            freq_obs = counts[r] / total_pair
            freq_ref = ref_counts[r] / total_ref
            
            # Score = -log(Obs/Ref) [cite: 32]
            if freq_obs > 0 and freq_ref > 0:
                ratio = freq_obs / freq_ref
                val = -math.log(ratio)
                scores.append(val)
            else:
                scores.append(10.0) # Max penalty [cite: 33]
                
        final_scores[pair] = scores
        
    return final_scores


def save_scores(scores_dict, output_dir="potentials"):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for pair, values in scores_dict.items():
        filename = os.path.join(output_dir, f"{pair}.txt")
        with open(filename, "w") as f:
            for v in values:
                f.write(f"{v:.4f}\n")
    print(f"Saved 10 scoring files to folder '{output_dir}'")

