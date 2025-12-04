
import math
import os
import glob
import argparse
import sys

from rna_loader import load_rna_structure
from rna_distance import get_all_distances
from rna_training import get_pair_name, get_all_distances, train_objective_function, save_scores

# Command-line input
parser = argparse.ArgumentParser(description="Train an RNA objective function from PDB or CIF files.")

parser.add_argument("-d", "--data", 
                    type=str, required=True, 
                    help="Path to the folder containing structure files (default is 'rna_data/pdb')")

parser.add_argument("-f", "--format", 
                    type=str, choices=["pdb", "cif"], default="pdb",
                    help="The file format to use for training (.pdb or cif). Default is pdb.")

parser.add_argument("-a", "--atom", 
                    type=str, default="C3'",
                    help="Atom type to calculate distances. Default is C3'.")

parser.add_argument("-m", "--mode", 
                    type=str, choices=["histogram", "kernel"], default="histogram",
                    help="Using histogram or Kernel density to plot the distances. Default is histogram.")

parser.add_argument("-b", "--bin_size", 
                    type=float, default=1.0,
                    help="Bin size for histogram")

parser.add_argument("-min_dist", "--min_dist", 
                    type=float, default=3.0,
                    help="Minimum distance between a pair of atom. Default is 3A.")

parser.add_argument("-max_dist", "--max_dist", 
                    type=float, default=20.0,
                    help="Maximum distance between a pair of atom. Default is 20A.")

parser.add_argument("-o", "--out_dir", 
                    type=str, default="potentials",
                    help="Output directory. Default is 'potentials'.")

args = parser.parse_args()


if __name__ == "__main__":
    # Validate input file directory
    if not os.path.isdir(args.data):
        print(f"Error: Input directory {args.data} does not exist.")
        sys.exit(1)

    # Find files based on user choice
    search_pattern = os.path.join(args.data, f"*.{args.format}")
    found_files = glob.glob(search_pattern)
    # also try uppercase (e.g., .PDB or .CIF) just in case
    if len(found_files) == 0:
        search_pattern_upper = os.path.join(args.data, f"*.{args.format.upper()}")
        found_files = glob.glob(search_pattern_upper)

    if len(found_files) == 0:
        print(f"No files in .{args.format} found in {args.data}")
        sys.exit(1)


    # Training
    print(f"Training on {len(found_files)} {args.format.upper()} files from {args.data}...")
    # scores = train_objective_function(found_files)
    scores = train_objective_function(found_files, 
                                      atom_type=args.atom, 
                                      mode=args.mode, 
                                      bin_size=args.bin_size, 
                                      max_dist=args.max_dist,
                                      min_dist=args.min_dist,
                                      bandwidth="scott")
    print(f"Saving scores to {args.out_dir}...")
    save_scores(scores, output_dir=args.out_dir)



