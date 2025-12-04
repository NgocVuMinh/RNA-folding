# rna_plot.py
# Script to plot RNA scoring curves from text files (correspoding to RNA base pair interactions).
# usage:
# KDE :  python rna_plot.py -i out_kernel -o plots_kde
# Histograms : python rna_plot.py -i out_hist -o plots_hist

import os
import argparse
import matplotlib.pyplot as plt

def plot_profiles(input_dir, output_dir):
    """
    Reads scoring profiles (.txt) and saves them as images (.png).
    """
    # 1. Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output folder: {output_dir}")

    # The 10 standard base pairs to look for
    pairs = ["AA", "AU", "AC", "AG", "UU", "UC", "UG", "CC", "CG", "GG"]
    
    # Store data to make a summary plot later
    all_data = {}

    print(f"Reading profiles from '{input_dir}'...")

    for pair in pairs:
        filename = f"{pair}.txt"
        filepath = os.path.join(input_dir, filename)
        
        if not os.path.exists(filepath):
            print(f"Skipping {pair} (file not found)")
            continue
            
        distances = []
        scores = []
        
        # 2. Read the file
        # Format is: Distance <tab> Score
        try:
            with open(filepath, "r") as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        distances.append(float(parts[0]))
                        scores.append(float(parts[1]))
            
            all_data[pair] = (distances, scores)

            # 3. Plot Individual Curve
            plt.figure(figsize=(8, 6))
            plt.plot(distances, scores, label=f"{pair} Interaction", color='blue', linewidth=2)
            
            # Formatting
            plt.axhline(0, color='black', linestyle='--', linewidth=1) # Zero energy line
            plt.title(f"Statistical Potential: {pair}")
            plt.xlabel("Distance (Å)")
            plt.ylabel("Pseudo-Energy Score")
            plt.grid(True, alpha=0.3)
            plt.legend()
            
            # Save individual plot
            out_name = os.path.join(output_dir, f"plot_{pair}.png")
            plt.savefig(out_name, dpi=300)
            plt.close() # Close memory to avoid warning
            print(f"Saved {out_name}")

        except Exception as e:
            print(f"Error processing {pair}: {e}")

    # 4. Create Summary Plot (All in one)
    if all_data:
        plt.figure(figsize=(12, 8))
        for pair, (dist, scr) in all_data.items():
            plt.plot(dist, scr, label=pair, linewidth=1.5, alpha=0.8)
            
        plt.axhline(0, color='black', linestyle='--', linewidth=1.5)
        plt.title("Comparison of All RNA Interactions")
        plt.xlabel("Distance (Å)")
        plt.ylabel("Pseudo-Energy Score")
        plt.ylim(-2, 5) # Zoom in on the interesting part (ignore huge penalties)
        plt.grid(True, alpha=0.3)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        
        summary_name = os.path.join(output_dir, "summary_all_pairs.png")
        plt.savefig(summary_name, dpi=300)
        print(f"Saved summary plot: {summary_name}")
    else:
        print("No data found to plot.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot RNA scoring curves from text files.")
    
    # Arguments
    parser.add_argument("-i", "--input", type=str, required=True, 
                        help="Input folder containing the .txt scoring files (e.g. 'potentials')")
    parser.add_argument("-o", "--output", type=str, default="plots", 
                        help="Output folder to save images (default: 'plots')")
    
    args = parser.parse_args()
    
    plot_profiles(args.input, args.output)