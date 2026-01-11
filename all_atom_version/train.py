# train.py (All-Atom Version)
import os
import glob
import argparse
import sys

# Import from local modules in this folder
from rna_training import train_objective_function_kernel, train_objective_function_histogram
from utils import save_scores, parse_bandwidth

# 1. Setup arguments
parser = argparse.ArgumentParser(description="Train an All-Atom RNA objective function.")

parser.add_argument("-d", "--data", type=str, required=True, 
                    help="Folder containing .pdb or .cif files")
parser.add_argument("-f", "--format", type=str, choices=["pdb", "cif"], default="pdb",
                    help="File format (default: pdb)")
# Atom argument removed/ignored since we use ALL atoms
parser.add_argument("-m", "--mode", type=str, choices=["histogram", "kernel"], default="histogram",
                    help="Scoring mode: 'histogram' or 'kernel' (default: histogram)")
parser.add_argument("-formula", "--formula", type=str, choices=["pmf", "tig"], default="tig",
                    help="Scoring function formula: 'pmf' or 'tig' (default: 'tig')")
parser.add_argument("-o", "--out_dir", type=str, default=None,
                    help="Output folder (default: profiles)")

# Advanced options
parser.add_argument("-b", "--bin_size", type=float, default=1.0, help="Histogram bin size")
parser.add_argument("-bw", "--bandwidth", type=parse_bandwidth, default="SJ", help="Bandwith for KDE")
parser.add_argument("-ktype", "--kernel_type", choices=["gaussian", "rectangular", "triangular", "epanechnikov", "biweight", "cosine", "optcosine", "gaussian"], default="gaussian", help="Kernel for KDE")
parser.add_argument("--min_dist", type=float, default=0.0, help="Min distance (A)")
parser.add_argument("--max_dist", type=float, default=20.0, help="Max distance (A)")
parser.add_argument("--plot", action="store_true", help="Enable plotting of distributions (WARNING: generates many files)")


args = parser.parse_args()

if __name__ == "__main__":
    # 2. Check input directory
    # Handle relative paths depending on where script is run
    data_path = os.path.abspath(args.data)
    
    if not os.path.isdir(data_path):
        print(f"Error: Directory '{data_path}' not found.")
        sys.exit(1)

    # 3. Find files
    pattern = os.path.join(data_path, f"*.{args.format}")
    found_files = glob.glob(pattern)
    
    # Try uppercase extension
    if not found_files:
        pattern_upper = os.path.join(data_path, f"*.{args.format.upper()}")
        found_files = glob.glob(pattern_upper)

    if not found_files:
        print(f"No .{args.format} files found in {data_path}")
        sys.exit(1)

    # Setting output directory:
    if args.out_dir is None:
        if args.mode=="histogram":
            out_dir = f"profiles/all_atom_hist_{args.formula}_bin{args.bin_size}"
        elif args.mode=="kernel":
            out_dir = f"profiles/all_atom_kde_{args.formula}_{args.kernel_type}_bw{args.bandwidth}"
    else:
        out_dir = args.out_dir

    # 4. Run training
    print(f"Starting All-Atom training on {len(found_files)} files...")
    print(f"Reading from: {data_path}")
    print(f"Output to:    {out_dir}")
    
    if args.mode=="histogram":
        try:
            scores = train_objective_function_histogram(
                found_files, 
                atom_type=None, # Ignored
                mode=args.mode, 
                bin_size=args.bin_size, 
                max_dist=args.max_dist,
                formula=args.formula,
                min_dist=args.min_dist,
                plot_dist=args.plot  # Pass the plot flag
            )
            
            save_scores(scores, output_dir=out_dir)
        
        except Exception as e:
            print(f"Error: {e}")
            import traceback
            traceback.print_exc()

    elif args.mode=="kernel":
        try:
            scores = train_objective_function_kernel(
                found_files, 
                atom_type=None, # Ignored
                mode=args.mode, 
                bandwidth=args.bandwidth, 
                kernel_type=args.kernel_type,
                max_dist=args.max_dist,
                formula=args.formula,
                min_dist=args.min_dist
            )
        
            save_scores(scores, output_dir=out_dir)
        
        except Exception as e:
            print(f"Error: {e}")
            import traceback
            traceback.print_exc()
