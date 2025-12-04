# utils.py
# This file contains utility functions for RNA structure processing

import os

def get_pair_name(res1, res2):
    """
    Standardizes pair names so order doesn't matter.
    Example: 'A' and 'U' -> 'AU'.
    """
    sorted_pair = sorted([res1, res2])
    return f"{sorted_pair[0]}{sorted_pair[1]}"

def save_scores(scores_dict, output_dir="outputs"):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    for pair, data_dict in scores_dict.items():
        pair_name = pair.replace(" ", "_").replace("/", "-")
        filename = os.path.join(output_dir, f"{pair_name}.txt")
        dists = data_dict['distances']
        scores = data_dict['scores']
        
        with open(filename, "w") as f:
            for d, s in zip(dists, scores):
                f.write(f"{d:.4f}\t{s:.4f}\n") # distance score
                
    print(f"Saved scoring profiles to folder '{output_dir}'")
