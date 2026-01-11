# utils.py
# All-atom version
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

def get_pair_name(res1, res2, atom1=None, atom2=None):
    """
    Standardizes pair names (alphabetical order).
    For all-atom version, can include atom names.
    e.g., 'U' and 'A' becomes 'AU'.
    If atom1/atom2 provided: 'U_C1prime' and 'A_N1' becomes 'AC1prime_AN1'
    """
    if atom1 is not None and atom2 is not None:
        pair1 = f"{res1}_{atom1}"
        pair2 = f"{res2}_{atom2}"
        sorted_pair = sorted([pair1, pair2])
        return f"{sorted_pair[0]}__{sorted_pair[1]}"
    else:
        sorted_pair = sorted([res1, res2])
        return f"{sorted_pair[0]}{sorted_pair[1]}"

def get_scoring_formula(formula):
    """
    Calculates the score given observed 
    and reference frequencies based on a formula PMF or TIG
    """
    if formula.lower() == "pmf":
        # potential mean of force -ln(obs / ref)
        return lambda obs, ref: -np.log(obs / ref)
    elif formula.lower() == "tig":
        # total information gain: -(obs - ref) / ref
        return lambda obs, ref: -(obs - ref) / ref
    else:
        raise ValueError(f"Unknown formula: {formula}")

def save_scores(scores_dict, output_dir="outputs"):
    """
    Saves scoring profiles to text files.
    Expects scores_dict[pair] = {'distances': [...], 'scores': [...]}
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    print(f"Saving {len(scores_dict)} profiles to '{output_dir}'...")
        
    for pair, data_dict in scores_dict.items():
        # Handle filename
        pair_name = pair.replace(" ", "_").replace("/", "-")
        filename = os.path.join(output_dir, f"{pair_name}.txt")
        
        # Extract data
        dists = data_dict['distances']
        scores = data_dict['scores']
        
        # Write to file
        with open(filename, "w") as f:
            for d, s in zip(dists, scores):
                f.write(f"{d:.4f}\t{s:.4f}\n")

def parse_bandwidth(x):
    if x in ["nrd0", "SJ", "nrd", "ucv", "bcv"]:
        return x
    try:
        return float(x)
    except ValueError:
        raise argparse.ArgumentTypeError("Bandwidth must be a float or 'nrd0', 'SJ', 'nrd', 'ucv', 'bcv'.")
    

def plot_distributions(pair_counts, ref_counts, bin_size, min_dist, max_dist, plot_dist_dir):
    """
    Plot individual histograms for each pair + the reference histogram.
    """
    if not os.path.exists(plot_dist_dir):
        os.makedirs(plot_dist_dir)

    num_bins = len(ref_counts)
    bin_centers = min_dist + np.arange(num_bins) * bin_size

    # XX (reference)
    plt.figure(figsize=(6, 4))
    plt.bar(bin_centers, ref_counts, width=bin_size*0.9)
    plt.xlabel("Distance (Å)")
    plt.ylabel("Count")
    plt.title("Reference (XX)")
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dist_dir, "XX.png"), dpi=300)
    plt.close()

    # Pairs
    for pair, counts in pair_counts.items():
        pname = "".join(pair.split("__")) if "__" in pair else "".join(pair)
        plt.figure(figsize=(6, 4))
        plt.bar(bin_centers, counts, width=bin_size * 0.9)
        plt.xlabel("Distance (Å)")
        plt.ylabel("Count")
        plt.title(f"{pair}")
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dist_dir, f"{pname}.png"), dpi=300)
        plt.close()
