# utils.py
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

def get_pair_name(res1, res2):
    """
    Standardizes pair names (alphabetical order).
    e.g., 'U' and 'A' becomes 'AU'.
    """
    sorted_pair = sorted([res1, res2])
    return f"{sorted_pair[0]}{sorted_pair[1]}"

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
    if x.lower() in ["scott", "silverman"]:
        return x.lower()
    try:
        return float(x)
    except ValueError:
        raise argparse.ArgumentTypeError("Bandwidth must be a float or 'scott'/'silverman'.")
    

def plot_distributions(pair_counts, ref_counts, bin_size, min_dist, max_dist, plot_dist_dir):
    """
    Plot individual histograms for each pair + the reference histogram.
    """
    if not os.path.exists(plot_dist_dir):
        os.makedirs(plot_dist_dir)

    num_bins = len(ref_counts)
    # bin_edges = np.arange(min_dist, max_dist + bin_size, bin_size)
    bin_centers = min_dist + np.arange(num_bins) * bin_size # bin_edges[:-1] + bin_size/2

    # XX
    plt.figure(figsize=(6, 4))
    plt.bar(bin_centers, ref_counts, width=bin_size*0.9)
    plt.xlabel("Distance (Å)")
    plt.ylabel("Count")
    plt.title("Reference (XX)")
    plt.tight_layout()
    # plt.show()
    plt.savefig(os.path.join(plot_dist_dir, "XX.png"), dpi=300)
    plt.close()

    # Pairs
    for pair, counts in pair_counts.items():
        pname="".join(pair)
        plt.figure(figsize=(6, 4))
        plt.bar(bin_centers, counts, width=bin_size * 0.9)
        plt.xlabel("Distance (Å)")
        plt.ylabel("Count")
        plt.title(f"{pair}")
        plt.tight_layout()
        #plt.show()
        plt.savefig(os.path.join(plot_dist_dir, f"{pname}.png"), dpi=300)
        plt.close()