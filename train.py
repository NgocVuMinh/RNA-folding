# main.py
import os
import glob
import argparse
import sys
from rna_training import train_objective_function
from utils import save_scores, parse_bandwidth

# 1. Setup arguments
parser = argparse.ArgumentParser(description="Train an RNA objective function.")

parser.add_argument("-d", "--data", type=str, required=True, 
                    help="Folder containing .pdb or .cif files")
parser.add_argument("-f", "--format", type=str, choices=["pdb", "cif"], default="pdb",
                    help="File format (default: pdb)")
parser.add_argument("-a", "--atom", type=str, default="C3'",
                    help="Atom name (default: C3')")
parser.add_argument("-m", "--mode", type=str, choices=["histogram", "kernel"], default="histogram",
                    help="Scoring mode: 'histogram' or 'kernel' (default: histogram)")
parser.add_argument("-o", "--out_dir", type=str, default="potentials",
                    help="Output folder (default: potentials)")
# Advanced options
parser.add_argument("-b", "--bin_size", type=float, default=1.0, help="Histogram bin size")
parser.add_argument("-bw", "--bandwidth", type=parse_bandwidth, default=1.0, help="Bandwith for KDE")
parser.add_argument("--min_dist", type=float, default=0.0, help="Min distance (A)")
parser.add_argument("--max_dist", type=float, default=20.0, help="Max distance (A)")

args = parser.parse_args()

if __name__ == "__main__":
    # 2. Check input directory
    if not os.path.isdir(args.data):
        print(f"Error: Directory '{args.data}' not found.")
        sys.exit(1)

    # 3. Find files
    pattern = os.path.join(args.data, f"*.{args.format}")
    found_files = glob.glob(pattern)
    
    # Try uppercase extension if lowercase found nothing
    if not found_files:
        pattern_upper = os.path.join(args.data, f"*.{args.format.upper()}")
        found_files = glob.glob(pattern_upper)

    if not found_files:
        print(f"No .{args.format} files found in {args.data}")
        sys.exit(1)

    # 4. Run training
    print(f"Starting training on {len(found_files)} files...")
    
    if args.mode=="histogram":
        try:
            scores = train_objective_function(
                found_files, 
                atom_type=args.atom, 
                mode=args.mode, 
                bin_size=args.bin_size, 
                max_dist=args.max_dist,
                min_dist=args.min_dist
            )
            
            # 5. Save results
            save_scores(scores, output_dir=args.out_dir)
        
        except Exception as e:
            print(f"Error: {e}")

    elif args.mode=="kernel":
        try:
            scores = train_objective_function(
                found_files, 
                atom_type=args.atom, 
                mode=args.mode, 
                bandwidth=args.bandwidth, 
                max_dist=args.max_dist,
                min_dist=args.min_dist
            )
        
            save_scores(scores, output_dir=args.out_dir)
        
        except Exception as e:
            print(f"Error: {e}")