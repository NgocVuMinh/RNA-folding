# rna_training.py
import numpy as np
import math
import os
from collections import defaultdict
# from sklearn.neighbors import KernelDensity
from scipy.stats import gaussian_kde
from rna_loader import load_rna_structure
from rna_distance import get_all_distances
from utils import get_pair_name, get_scoring_formula, plot_distributions

def train_objective_function_histogram(structure_files, 
                             atom_type=None,  # Ignored in all-atom
                             mode="histogram", 
                             bin_size=1.0,
                             min_dist=0.0, 
                             formula="tig",
                             seq_sep=3, 
                             max_dist=20.0,
                             plot_dist=False, # Default to False to avoid generating hundreds of plots
                             plot_dist_dir="distance_distributions_all_atom"):
    """
    Trains the statistical potential using histogram
    
    Returns:
        dict: { 'A_N1__U_O2': {'distances': [x...], 'scores': [y...]} ... }
    """
    valid_bases = ['A', 'U', 'C', 'G']
    
    # --- 1. Initialization ---
    num_bins = int(max_dist)
    
    # Use defaultdict for dynamic pair discovery (e.g. A_N1__U_O2)
    # Mapping pair_name -> dict of bin_index -> count
    pair_counts = defaultdict(lambda: np.zeros(num_bins, dtype=int))
    
    # Reference (XX) counts
    ref_counts = np.zeros(num_bins, dtype=int)
    
    print(f"Training on {len(structure_files)} structures (All-Atom mode, bin_size={bin_size}, max_dist={max_dist})")

    # --- 2. Data collection ---
    total_interactions = 0
    
    for filepath in structure_files:
        structure = load_rna_structure(filepath)
        if structure is None:
            continue
        
        interactions = get_all_distances(structure, 
                                         atom_name=None, # Use all atoms
                                         max_distance=max_dist,
                                         seq_sep=seq_sep)
        
        for item in interactions:
            if item["Type"] != "Intrachain":
                continue
            
            r1, r2 = item['Res1'], item['Res2']
            if r1 not in valid_bases or r2 not in valid_bases:
                continue
            
            # Use atom-specific naming
            pname = get_pair_name(r1, r2, item['AtomA'], item['AtomB'])
            
            dist = item['Distance']
            
            bin_idx = int(math.floor(dist))
            if 0 <= bin_idx < len(ref_counts):
                pair_counts[pname][bin_idx] += 1
                ref_counts[bin_idx] += 1
                total_interactions += 1

    if plot_dist:
        # Note: This might generate A LOT of files in all-atom mode
        print(f"Plotting distributions for {len(pair_counts)} pairs (this may take a while)...")
        plot_distributions(pair_counts, ref_counts, bin_size, min_dist, max_dist, plot_dist_dir)
        print(f"Distance distributions saved to {plot_dist_dir}.")

    if total_interactions == 0:
        print("Warning: No valid interactions found.")
        return {}

    print(f"Total interactions processed: {total_interactions}")
    print(f"Unique atom pairs found: {len(pair_counts)}")
    print("Computing scores...")
    
    final_scores = {}

    # --- 3. Scoring ---
    calculate_score = get_scoring_formula(formula)
    
    total_ref = ref_counts.sum()
    freq_ref = ref_counts / total_ref
    
    # Iterate over whatever pairs we found
    for pair, counts in pair_counts.items():
        scores = []
        total_pair = counts.sum()
        
        if total_pair == 0:
            final_scores[pair] = {'distances': list(range(num_bins)),
                                  'scores': [10.0] * num_bins}
            continue
        
        freq_obs = counts / total_pair
        
        for i in range(len(counts)):
            if freq_obs[i] > 0 and freq_ref[i] > 0:
                score = calculate_score(freq_obs[i], freq_ref[i])
                scores.append(score)
            else:
                scores.append(10.0)
        
        final_scores[pair] = {
            'distances': list(range(num_bins)),
            'scores': scores
        }
            
    return final_scores



def train_objective_function_kernel(structure_files, 
                             atom_type=None, 
                             mode="kernel", 
                             min_dist=0.0, 
                             max_dist=20.0,
                             formula="tig",
                             seq_sep=3,
                             kernel_type="gaussian",
                             bandwidth="SJ"):
    """
    Trains the statistical potential using Kernel Density Estimation (All-Atom)
    Uses scipy.stats.gaussian_kde (simpler alternative to R).
    """
    valid_bases = ['A', 'U', 'C', 'G']
    
    # --- 1. Initialization ---
    # Store raw distances for KDE
    pair_counts = defaultdict(list)
    ref_counts = []
    
    # Resolve bandwidth for Scipy
    # Scipy accepts 'scott', 'silverman' or a scalar
    scipy_bw = 'scott' # Default
    if isinstance(bandwidth, (int, float)):
        scipy_bw = bandwidth
    elif bandwidth in ['nrd0', 'nrd', 'ucv', 'bcv', 'SJ']:
        # Map R-specific selectors to 'scott' (closest approximation)
        print(f"Note: Bandwidth selector '{bandwidth}' is R-specific. Using 'scott' rule for scipy.")
        scipy_bw = 'scott'
    elif bandwidth in ['scott', 'silverman']:
        scipy_bw = bandwidth

    if kernel_type != 'gaussian':
         print(f"Note: Kernel '{kernel_type}' is not supported by scipy (gaussian only). Using gaussian.")

    print(f"Training on {len(structure_files)} structures (All-Atom mode, kernel=gaussian, bw={scipy_bw})")

    # --- 2. Data collection ---
    total_interactions = 0
    
    for filepath in structure_files:
        structure = load_rna_structure(filepath)
        if structure is None:
            continue
        
        interactions = get_all_distances(structure, 
                                         atom_name=None,
                                         max_distance=max_dist,
                                         seq_sep=seq_sep)
        
        for item in interactions:
            if item["Type"] != "Intrachain":
                continue
            
            r1, r2 = item['Res1'], item['Res2']
            if r1 not in valid_bases or r2 not in valid_bases:
                continue
            
            dist = item['Distance']
            pname = get_pair_name(r1, r2, item['AtomA'], item['AtomB'])
            
            pair_counts[pname].append(dist)
            ref_counts.append(dist)
            total_interactions += 1

    if total_interactions == 0:
        print("Warning: No valid interactions found.")
        return {}

    print(f"Total interactions processed: {total_interactions}")
    print(f"Unique atom pairs found: {len(pair_counts)}")
    print("Computing scores...")
    
    final_scores = {}

    # --- 3. Scoring ---
    calculate_score = get_scoring_formula(formula)
    
    if not ref_counts: 
        print("Warning: reference distance list is empty")
        return {}
        
    x_grid = np.linspace(min_dist, max_dist, 200)
    
    # Reference KDE (XX)
    try:
        ref_kde = gaussian_kde(ref_counts, bw_method=scipy_bw)
        ref_pdf = ref_kde.evaluate(x_grid)
        ref_pdf = np.maximum(ref_pdf, 1e-6) # Avoid 0
    except Exception as e:
        print(f"Error computing reference KDE: {e}")
        return {}
    
    for pair, dist_values in pair_counts.items():
        if len(dist_values) < 2:
            final_scores[pair] = {'distances': x_grid.tolist(),
                                    'scores': [10.0] * 200}
            continue
        
        # Pair KDE (XY)
        try:
            pair_kde = gaussian_kde(dist_values, bw_method=scipy_bw)
            obs_pdf = pair_kde.evaluate(x_grid)
            obs_pdf = np.maximum(obs_pdf, 1e-6)
            
            scores = calculate_score(obs_pdf, ref_pdf)
            # Forcing max penalty for very short distances (clashes)
            scores[x_grid < 2.0] = 10.0
            
            final_scores[pair] = {'distances': x_grid.tolist(),
                                  'scores': scores.tolist()}
        except Exception as e:
            print(f"Error computing KDE for {pair}: {e}")
            final_scores[pair] = {'distances': x_grid.tolist(),
                                  'scores': [10.0] * 200}
            
    return final_scores