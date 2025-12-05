import os

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