# RNA Training Script
# This script trains an RNA objective function based on provided PDB or CIF files.


import numpy as np
import math
import os
import glob
import argparse
import sys
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
    Train distance-based statistical potentials for RNA structure evaluation.
    
    Output: dict, Dictionary mapping base pairs to score arrays/functions
    """
    
    valid_bases = ['A', 'U', 'C', 'G']
    
    
    if mode == "histogram": # discrete bins
        num_bins = int((max_dist - min_dist) / bin_size) + 1
        bin_edges = np.linspace(min_dist, max_dist, num_bins)
        
        pair_counts = {}
        ref_counts = np.zeros(num_bins - 1, dtype=float)
        
        for i in range(len(valid_bases)):
            for j in range(i, len(valid_bases)):
                pname = get_pair_name(valid_bases[i], valid_bases[j])
                pair_counts[pname] = np.zeros(num_bins - 1, dtype=float)
        
        print(f"Constructing histograms: {num_bins-1} bins from {min_dist}-{max_dist}A...")
    
    elif mode == "kernel": # continious raw distances
        pair_raw = {}
        ref_raw = []
        
        for i in range(len(valid_bases)):
            for j in range(i, len(valid_bases)):
                pname = get_pair_name(valid_bases[i], valid_bases[j])
                pair_raw[pname] = []
        
        print(f"Constructing Kernel density...")
    
    print(f"Training on {len(structure_files)} structure files...")
    
    total_interactions = 0
    
    for idx, filepath in enumerate(structure_files, 1):
        
        structure = load_rna_structure(filepath)
        if structure is None:
            print(f"Warning: Structure {filepath} is empty.")
            continue
        
        interactions = get_all_distances(structure,
                                         atom_name=atom_type,
                                         min_distance=min_dist,
                                         max_distance=max_dist)
        
        file_count = 0
        for item in interactions:
            if item["Type"] != "Intrachain":
                continue
            
            r1 = item['Res1']
            r2 = item['Res2']
            dist = item['Distance']
            
            if r1 not in valid_bases or r2 not in valid_bases:
                continue
            
            pname = get_pair_name(r1, r2)
            
            if mode == "histogram":
                bin_idx = np.searchsorted(bin_edges, dist) - 1
                if 0 <= bin_idx < len(ref_counts):
                    pair_counts[pname][bin_idx] += 1
                    ref_counts[bin_idx] += 1
                    file_count += 1
            
            elif mode == "kernel":
                pair_raw[pname].append(dist)
                ref_raw.append(dist)
                file_count += 1
        
        total_interactions += file_count
    
    if total_interactions == 0:
        raise ValueError("No valid interactions found in all structures.")
    
    print(f"\nTotal number of interactions: {total_interactions}")
    
    print("\nComputing scores...")
    final_scores = {}
    

    # ----- Histogram -----
    if mode == "histogram":
        total_ref = ref_counts.sum()
        if total_ref == 0:
            raise ValueError("Distribution is 0.")
        
        freq_ref = ref_counts / total_ref
        
        for pair, counts in pair_counts.items():
            scores = []
            total_pair = counts.sum()
            
            if total_pair == 0:
                final_scores[pair] = np.full(len(counts), 10.0).tolist() # penalty=10
                continue
            
            freq_obs = counts / total_pair
            
            # Score = -log(Obs/Ref)
            for i in range(len(counts)):
                if freq_obs[i] > 0 and freq_ref[i] > 0:
                    ratio = freq_obs[i] / freq_ref[i]
                    score = -math.log(ratio)
                    scores.append(score)
                    # score = -math.log(freq_obs[i] / freq_ref[i])
                    # score = min(score, 10.0)
                    # score = max(score, -10.0)
                else:
                    scores.append(10.0) # penalty=10 for no observation in a bin
            
            final_scores[pair] = scores

    
    # ----- Kernel density -----
    elif mode == "kernel":
        if len(ref_raw) == 0:
            raise ValueError("No distances collected.")
        
        # reference = standard gaussian kde
        ref_kde = gaussian_kde(ref_raw, bw_method="scott")
        eval_points = np.linspace(min_dist, max_dist, 100)
        ref_pdf = ref_kde(eval_points)
        ref_pdf = np.maximum(ref_pdf, 1e-6) # prevent zeros
        
        for pair, dist_values in pair_raw.items():
            if len(dist_values) == 0:
                scores = np.full(len(eval_points), 10.0).tolist()
                # final_scores[pair] = {
                #     'distances': eval_points.tolist(),
                #     'scores': np.full(len(eval_points), 10.0).tolist()
                # }
                print(f"Warning: No {pair} contacts found?")
                continue
            
            print(f"Fitting KDE for {len(dist_values)} {pair}...")
            pair_kde = gaussian_kde(dist_values, bw_method=bandwidth)
            obs_pdf = pair_kde(eval_points)
            obs_pdf = np.maximum(obs_pdf, 1e-6)
            
            scores = -np.log(obs_pdf / ref_pdf)
            # scores = np.clip(scores, -10.0, 10.0) # soft clamping
            
            final_scores[pair] = {'distances': eval_points.tolist(),
                                  'scores': scores.tolist()}
    
    print(f"\nTraining {mode} completed.")
    
    return final_scores


