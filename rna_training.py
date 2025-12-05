import numpy as np
import math
import os
from scipy.stats import gaussian_kde
from rna_loader import load_rna_structure
from rna_distance import get_all_distances
from utils import get_pair_name

def train_objective_function(structure_files, 
                             atom_type="C3'", 
                             mode="histogram", 
                             bin_size=1.0,
                             min_dist=3.0, 
                             max_dist=20.0,
                             bandwidth="scott"):
    """
    Trains the statistical potential using PDB/CIF files.
    
    Returns:
        dict: { 'AU': {'distances': [x...], 'scores': [y...]} ... }
    """
    valid_bases = ['A', 'U', 'C', 'G']
    
    # --- 1. Initialization ---
    if mode == "histogram":
        # Create bins (e.g., 3.0, 4.0, ... 20.0)
        num_bins = int((max_dist - min_dist) / bin_size) + 1
        bin_edges = np.linspace(min_dist, max_dist, num_bins)
        # Calculate centers for plotting (x-axis)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        pair_counts = {}
        # Reference (XX) counts
        ref_counts = np.zeros(num_bins - 1, dtype=float)
        
        # Initialize dictionary for all 10 pairs
        for i in range(len(valid_bases)):
            for j in range(i, len(valid_bases)):
                pname = get_pair_name(valid_bases[i], valid_bases[j])
                pair_counts[pname] = np.zeros(num_bins - 1, dtype=float)
                
    elif mode == "kernel":
        # Store raw distances for KDE
        pair_raw = {}
        ref_raw = []
        for i in range(len(valid_bases)):
            for j in range(i, len(valid_bases)):
                pname = get_pair_name(valid_bases[i], valid_bases[j])
                pair_raw[pname] = []

    print(f"Training on {len(structure_files)} files (Mode: {mode})...")

    # --- 2. Data Collection ---
    total_interactions = 0
    
    for filepath in structure_files:
        structure = load_rna_structure(filepath)
        if structure is None:
            continue
        
        interactions = get_all_distances(structure, atom_name=atom_type,
                                         min_distance=min_dist, max_distance=max_dist)
        
        for item in interactions:
            # Training uses Intrachain only
            if item["Type"] != "Intrachain":
                continue
            
            r1, r2 = item['Res1'], item['Res2']
            if r1 not in valid_bases or r2 not in valid_bases:
                continue
            
            dist = item['Distance']
            pname = get_pair_name(r1, r2)
            
            if mode == "histogram":
                # Find which bin this distance falls into
                bin_idx = np.searchsorted(bin_edges, dist) - 1
                if 0 <= bin_idx < len(ref_counts):
                    pair_counts[pname][bin_idx] += 1
                    ref_counts[bin_idx] += 1
                    total_interactions += 1
            
            elif mode == "kernel":
                pair_raw[pname].append(dist)
                ref_raw.append(dist)
                total_interactions += 1

    if total_interactions == 0:
        print("Warning: No valid interactions found.")
        return {}

    print(f"Total interactions processed: {total_interactions}")
    print("Computing scores...")
    
    final_scores = {}

    # --- 3. Scoring (Inverse Boltzmann) ---
    # Score = -ln( P_obs / P_ref )
    
    if mode == "histogram":
        total_ref = ref_counts.sum()
        freq_ref = ref_counts / total_ref
        
        for pair, counts in pair_counts.items():
            scores = []
            total_pair = counts.sum()
            
            if total_pair == 0:
                # No data for this pair -> Max penalty
                final_scores[pair] = {'distances': bin_centers.tolist(), 
                                      'scores': [10.0] * len(bin_centers)}
                continue
            
            freq_obs = counts / total_pair
            
            for i in range(len(counts)):
                # Avoid division by zero or log(0)
                if freq_obs[i] > 0 and freq_ref[i] > 0:
                    ratio = freq_obs[i] / freq_ref[i]
                    score = -math.log(ratio)
                    scores.append(score)
                else:
                    scores.append(10.0) # Penalty for unobserved bins
            
            # FIX: Return dictionary matching KDE format
            final_scores[pair] = {
                'distances': bin_centers.tolist(),
                'scores': scores
            }

    elif mode == "kernel":
        if not ref_raw: 
            return {}
            
        # Reference KDE (XX)
        ref_kde = gaussian_kde(ref_raw, bw_method="scott")
        eval_points = np.linspace(min_dist, max_dist, 100)
        ref_pdf = ref_kde(eval_points)
        ref_pdf = np.maximum(ref_pdf, 1e-10) # Avoid 0
        
        for pair, dist_values in pair_raw.items():
            if len(dist_values) < 2: # Need data points for KDE
                final_scores[pair] = {'distances': eval_points.tolist(),
                                      'scores': [10.0] * len(eval_points)}
                continue
            
            # Pair KDE (XY)
            pair_kde = gaussian_kde(dist_values, bw_method=bandwidth)
            obs_pdf = pair_kde(eval_points)
            obs_pdf = np.maximum(obs_pdf, 1e-10)
            
            # Vectorized scoring
            scores = -np.log(obs_pdf / ref_pdf)
            
            final_scores[pair] = {
                'distances': eval_points.tolist(),
                'scores': scores.tolist()
            }
            
    return final_scores