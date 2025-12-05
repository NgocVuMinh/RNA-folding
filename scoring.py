
import numpy as np
import os
import argparse
from scipy.interpolate import interp1d
from rna_loader import load_rna_structure
from rna_distance import is_strict_pure_rna, get_all_distances


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


def score_structure(structure_file, profile_dir, atom_name="C3'", 
                   max_distance=20.0, seq_sep=4):
    """Score an RNA structure."""

    valid_bases = ['A', 'U', 'C', 'G']

    # Load inputs
    profiles = load_profiles(profile_dir)
    structure = load_rna_structure(structure_file)
    interpolators = linear_interpol(profiles)
    
    # Compute distances
    interactions = get_all_distances(structure, atom_name, max_distance, seq_sep=3)
    
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
        if pair in interpolators:
            score = float(interpolators[pair](dist))
        else:
            print(f"Warning: No potential for {pair}, max penalty=10")
            score = 10.0
        
        total_energy += score
    
    results = {'structure': structure_file,
               'total_energy': total_energy}
    
    return results

    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Score an RNA structure using trained profiles')
    parser.add_argument( "--structure_dir", type=str, required=True,
                        help="RNA structure file (PDB or CIF)")
    parser.add_argument("--profile_dir", type=str, required=True,
                        help="Directory containing potential files")
    parser.add_argument("--atom_name", type=str, default="C3'",
                        help="Atom name (default: C3')")
    parser.add_argument("--max_distance", type=float, default=20.0, 
                        help="Max distance (A)")
    
    args = parser.parse_args()
    
    results = score_structure(args.structure_dir,
                              args.profile_dir,
                              atom_name=args.atom_name,
                              max_distance=args.max_distance)
    
    if results:
        print(f"Structure: {results['structure']}")
        print(f"Estimated Gibbs free energy: {results['total_energy']:.3f}")