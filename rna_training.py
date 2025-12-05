# rna_training.py
import numpy as np
import math
import os
# from scipy.stats import gaussian_kde
from sklearn.neighbors import KernelDensity
from rna_loader import load_rna_structure
from rna_distance import get_all_distances
from utils import get_pair_name, plot_distributions

# to use stats.density from R 
from rpy2.robjects.packages import importr
from rpy2.robjects import vectors
stats = importr("stats")

def train_objective_function_histogram(structure_files, 
                             atom_type="C3'", 
                             mode="histogram", 
                             bin_size=1.0,
                             min_dist=0.0, 
                             seq_sep=3, # only consider residues separated by at least 3 positions on the sequence 
                             max_dist=20.0,
                             plot_dist=True,
                             plot_dist_dir="distance_distributions"):
    """
    Trains the statistical potential using histogram
    
    Returns:
        dict: { 'AU': {'distances': [x...], 'scores': [y...]} ... }
    """
    valid_bases = ['A', 'U', 'C', 'G']
    
    # --- 1. Initialization ---
    # Create bins (e.g., 3.0, 4.0, ... 20.0)
    num_bins = int((max_dist - min_dist) / bin_size) # int(max_dist)
    # bin_edges = np.linspace(min_dist, max_dist, num_bins)
    # Calculate centers for plotting (x-axis)
    # bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    pair_counts = {}
    # Reference (XX) counts
    ref_counts = np.zeros(num_bins, dtype=int)
    
    # Initialize dictionary for all 10 pairs
    for i in range(len(valid_bases)):
        for j in range(i, len(valid_bases)):
            pname = get_pair_name(valid_bases[i], valid_bases[j])
            pair_counts[pname] = np.zeros(num_bins, dtype=int)
                

    print(f"Training on {len(structure_files)} structures (atom: {atom_type}, mode: {mode}, bin_size={bin_size}, min_dist={min_dist}, max_dist={max_dist})")

    # --- 2. Data collection ---
    total_interactions = 0
    
    for filepath in structure_files:
        structure = load_rna_structure(filepath)
        if structure is None:
            continue
        
        interactions = get_all_distances(structure, 
                                         atom_name=atom_type,
                                         # min_distance=min_dist,
                                         max_distance=max_dist,
                                         seq_sep=seq_sep)
        
        for item in interactions:
            # Training uses Intrachain only
            if item["Type"] != "Intrachain":
                continue
            
            r1, r2 = item['Res1'], item['Res2']
            if r1 not in valid_bases or r2 not in valid_bases:
                continue
            
            dist = item['Distance']
            pname = get_pair_name(r1, r2)
            
            # Find which bin this distance falls into
            bin_idx = int(math.floor(dist)) #np.searchsorted(bin_edges, dist) - 1
            if 0 <= bin_idx < len(ref_counts):
                pair_counts[pname][bin_idx] += 1
                ref_counts[bin_idx] += 1
                total_interactions += 1

    if plot_dist:
        plot_distributions(pair_counts, ref_counts, bin_size, min_dist, max_dist, plot_dist_dir)
        print(f"Distance distributions saved to {plot_dist_dir}.")

    if total_interactions == 0:
        print("Warning: No valid interactions found.")
        return {}

    print(f"Total interactions processed: {total_interactions}")
    print("Computing scores...")
    
    final_scores = {}

    # --- 3. Scoring (Inverse Boltzmann) ---
    # Score = -ln( P_obs / P_ref )
    
    total_ref = ref_counts.sum()
    freq_ref = ref_counts / total_ref
    
    for pair, counts in pair_counts.items():
        scores = []
        total_pair = counts.sum()
        
        if total_pair == 0:
            # No data for this pair -> Max penalty
            final_scores[pair] = {'distances': list(range(num_bins)), 
                                    'scores': [10.0] * num_bins}
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
            'distances': list(range(num_bins)),
            'scores': scores
        }
            
    return final_scores



def train_objective_function_kernel(structure_files, 
                             atom_type="C3'", 
                             mode="kernel", 
                             min_dist=0.0, max_dist=20.0,
                             seq_sep=3, # only consider residues separated by at least 3 positions on the sequence 
                             kernel_type="gaussian",
                             bandwidth="SJ"):
    """
    Trains the statistical potential using PDB/CIF files.
    
    Returns:
        dict: { 'AU': {'distances': [x...], 'scores': [y...]} ... }
    """
    valid_bases = ['A', 'U', 'C', 'G']
    
    # --- 1. Initialization ---
    # Store raw distances for KDE
    pair_counts = {}
    ref_counts = []
    for i in range(len(valid_bases)):
        for j in range(i, len(valid_bases)):
            pname = get_pair_name(valid_bases[i], valid_bases[j])
            pair_counts[pname] = []

    print(f"Training on {len(structure_files)} structures (atom: {atom_type}, mode: {mode}, kernel_type={kernel_type}, bandwidth={bandwidth})")

    # --- 2. Data collection ---
    total_interactions = 0
    
    for filepath in structure_files:
        structure = load_rna_structure(filepath)
        if structure is None:
            continue
        
        interactions = get_all_distances(structure, 
                                         atom_name=atom_type,
                                         # min_distance=min_dist,
                                         max_distance=max_dist,
                                         seq_sep=seq_sep)
        
        for item in interactions:
            # Training uses Intrachain only
            if item["Type"] != "Intrachain":
                continue
            
            r1, r2 = item['Res1'], item['Res2']
            if r1 not in valid_bases or r2 not in valid_bases:
                continue
            
            dist = item['Distance']
            pname = get_pair_name(r1, r2)
            
            pair_counts[pname].append(dist)
            ref_counts.append(dist)
            total_interactions += 1

    if total_interactions == 0:
        print("Warning: No valid interactions found.")
        return {}

    print(f"Total interactions processed: {total_interactions}")
    print("Computing scores...")
    
    final_scores = {}

    # --- 3. Scoring (Inverse Boltzmann) ---
    # Score = -ln( P_obs / P_ref )
    
    if not ref_counts: 
        print("Warning: reference distance list is empty")
        return {}
        
    # Reference KDE (XX)
    ref_counts = vectors.FloatVector(ref_counts)
    ref_kde = stats.density(ref_counts, 
                            adjust=1,
                            bw=bandwidth, kernel=kernel_type,
                            from_=min_dist, to=max_dist, n=200)
    ref_pdf = np.array(ref_kde.rx2("y"))
    ref_pdf = np.maximum(ref_pdf, 1e-6) # Avoid 0
    
    for pair, dist_values in pair_counts.items():
        if len(dist_values) < 2: # Need data points for KDE
            final_scores[pair] = {'distances': np.linspace(min_dist, max_dist, 200).tolist(),
                                    'scores': [10.0] * 200}
            continue
        
        # Pair KDE (XY)
        dist_values = vectors.FloatVector(dist_values)
        pair_kde = stats.density(dist_values, 
                                 adjust=1,
                                 bw=bandwidth, kernel=kernel_type,
                                 from_=min_dist, to=max_dist, n=200)
        obs_pdf = np.array(pair_kde.rx2("y"))
        obs_pdf = np.maximum(obs_pdf, 1e-6)
        
        scores = -np.log(obs_pdf / ref_pdf)
        
        final_scores[pair] = {'distances': np.linspace(min_dist, max_dist, 200).tolist(),
                              'scores': scores.tolist()}
            
    return final_scores