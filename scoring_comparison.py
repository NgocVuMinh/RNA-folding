import numpy as np
import pandas as pd
import os, sys, glob
import argparse
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
import itertools

from rna_loader import load_rna_structure
from rna_distance import get_all_distances


def load_profiles(profile_dir):
    """Load trained profiles from text files."""
    
    pairs = ["AA", "AU", "AC", "AG", "UU", "CU", "GU", "CC", "CG", "GG"]
    profiles = {}
    
    print(f"Reading profiles from '{profile_dir}'...")
    
    for pair in pairs:
        filename = f"{pair}.txt"
        filepath = os.path.join(profile_dir, filename)
        
        if not os.path.exists(filepath):
            print(f"Skipping {pair} (file not found)")
            continue
        
        distances = []
        scores = []
        
        # Read the file (format: distance<tab>score)
        with open(filepath, "r") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    distances.append(float(parts[0]))
                    scores.append(float(parts[1]))
        
        # pair string to tuple (e.g., "AU" -> ('A', 'U'))
        pair_tuple = tuple(sorted([pair[0], pair[1]]))
        
        profiles[pair_tuple] = {'distances': distances,
                                'scores': scores}
    
    return profiles


def linear_interpol(profiles):
    """Linear interpolation"""
    interpolators = {}
    
    for pair, data in profiles.items():
        distances = np.array(data['distances'])
        scores = np.array(data['scores'])
        
        interpolators[pair] = interp1d(distances,
            scores,
            kind='linear',
            bounds_error=False,
            fill_value=(scores[0], scores[-1]))
    
    return interpolators


def score_structure(structure, 
                    interpolator, 
                    atom_name="C3'", 
                    max_distance=20.0, seq_sep=3):
    """Score an RNA structure."""

    valid_bases = {'A', 'U', 'C', 'G'}

    # Load inputs
    # profile = load_profiles(profile_dir)
    # structure = load_rna_structure(structure_file)
    # interpolators = linear_interpol(profile)
    
    # Compute distances
    interactions = get_all_distances(structure, atom_name, max_distance, seq_sep=seq_sep)
    
    if len(interactions) == 0:
        print("Warning: No valid interactions found.")
        return None
    
    # Score all interactions
    total_energy = 0.0
    
    for item in interactions:
        if item["Type"] != "Intrachain":
                continue
        
        r1, r2 = item['Res1'], item['Res2']
        if r1 not in valid_bases or r2 not in valid_bases:
                continue
        
        dist = item['Distance']
        pair = tuple(sorted([r1, r2])) # pname = get_pair_name(r1, r2)
        
        # Interpolated score
        if pair in interpolator:
            score = float(interpolator[pair](dist))
        else:
            # print(f"Warning: No potential for {pair}, max penalty=10")
            score = 10.0
        
        total_energy += score
    
    # results = {'structure': structure_file,
    #            'total_energy': total_energy}
    
    return total_energy

    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Score an RNA structure using trained profiles')
    parser.add_argument( "--structure_dir", type=str, required=True,
                        help="RNA structure file (PDB or CIF)")
    parser.add_argument("--profile_dir", type=str, required=True,
                        nargs="+", # allowing multple profiles to be loaded
                        help="Directory containing potential files")
    parser.add_argument( "--out_dir", type=str, default="scoring_comparison",
                        help="Output path to store scoring results of multiple profiles on multiple structures.")
    parser.add_argument("--atom_name", type=str, default="C3'",
                        help="Atom name (default: C3')")
    parser.add_argument("--max_distance", type=float, default=20.0, 
                        help="Max distance (A)")
    
    args = parser.parse_args()


    # --- Loading profiles and interpolating
    interpolators = {}
    for prf in args.profile_dir:
        profile = load_profiles(prf)
        interpolators[prf] = linear_interpol(profile)
    print(f"Number of profiles loaded: {len(interpolators)}.")


    # --- Reading RNA structures
    if not os.path.isdir(args.structure_dir):
        print(f"Error: Structure directory '{args.structure_dir}' not found.")
        sys.exit(1)

    print(f"Loading RNA structures from {args.structure_dir}...")
    rna_files = glob.glob(os.path.join(args.structure_dir, f"*.pdb")) + glob.glob(os.path.join(args.structure_dir, f"*.cif"))

    if len(rna_files) == 0:
        print(f"Error: Structure directory '{args.structure_dir}' is empty.")
        sys.exit(1)


    # --- Scoring
    energy_results = []
    for rna_file in rna_files:
        structure = load_rna_structure(rna_file)
        rna_name = os.path.basename(rna_file)
        
        res = {}
        for profile_dir, interpol in interpolators.items():
            profile_name = os.path.basename(os.path.normpath(profile_dir))
            energy = score_structure(structure,
                                   interpol,
                                   atom_name=args.atom_name,
                                   max_distance=args.max_distance)
            if energy is None:
                print(f"Energy is None for {rna_name}")
                continue
            else:
                energy_results.append({"rna": rna_name,
                                   "profile": profile_name,
                                   "energy": energy})
        
    # --- Saving scoring outputs
    os.makedirs(args.out_dir, exist_ok=True)
    with open(os.path.join(args.out_dir, "energies.txt"), "w") as out:
        out.write("rna\tprofile\tenergy\n")
        for row in energy_results:
            out.write(f"{row['rna']}\t{row['profile']}\t{row['energy']}\n")

    # --- Ploting the correlation between different profiles (if there are more than 1 profile)
    if len(args.profile_dir) > 1:
        df = pd.DataFrame(energy_results)
        df = df[df["energy"] > -10] # we remove extreme negative scores (potentially of very large structures)
        energy_matrix = df.pivot(index="rna", columns="profile", values="energy")
        

        profiles = energy_matrix.columns.tolist()

        for p1, p2 in itertools.combinations(profiles, 2):
            x = energy_matrix[p1]
            y = energy_matrix[p2]

            # drop RNAs missing either score
            mask = x.notna() & y.notna()
            x = x[mask]
            y = y[mask]

            if len(x) < 2:
                continue  # not enough points to correlate

            corr = x.corr(y, method="pearson")

            plt.figure(figsize=(4, 4))
            plt.scatter(x, y, alpha=0.7, s=8, color="crimson")
            plt.xlabel(p1)
            plt.ylabel(p2)
            plt.title(f"Pearson r={corr:.4f}")
            plt.tight_layout()

            plt.savefig(os.path.join(args.out_dir, f"{p1}_VS_{p2}.png"), dpi=300)
            plt.close()