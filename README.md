# RNA-folding

**M2 GENIOMHE - Univ Evry Paris-Saclay (2025-2026)**

**Our project for the "Bioinformatics of RNA and non-coding world" course**

- **Minh Ngoc VU**
- **Erine Benoist**

## Overview

For RNA structures, a native folding is the one with the lowest Gibbs free energy. The purpose of the project is to develop an objective function to estimate the energy and thus evaluating how close a predicted fold is to the optimal / native fold. 

The scoring function is a statistical potential derived from known experimental data. It calculates pseudo-energies based on the probability of observing specific atomic distances between residue pairs.

## Dataset

We obtained 100 RNA structures from *Homo sapiens* from the [RCSB PDB REST API](https://www.rcsb.org/docs/programmatic-access/web-apis-overview). Each structure was obtained in both PDB and CIF formats.

## Loading PDB and mmCIF files

The project utilizes the Biopython library (Bio.PDB) to parse structural data. The loading module (rna_loader.py) is designed to be format-agnostic:

- Format Detection: The script automatically detects the file extension (.pdb or .cif) to instantiate the correct parser (PDBParser or MMCIFParser).

- Model Selection: For files containing multiple models (common in NMR structures), the loader automatically selects the first model (Model 0) to ensure consistency across the dataset.

## Cleaning the Chains

To generate a high-quality statistical potential, we implemented a strict filtering mechanism in rna_distance.py. Experimental PDB files often contain "noise" such as protein chains, DNA-RNA hybrids, or ambiguous IUPAC codes.

Our is_strict_pure_rna function applies the following logic:

- Standard Bases Only: It scans every chain to ensure it is composed exclusively of standard RNA bases (Adenine, Uracil, Cytosine, Guanine).

- Whole-Chain Rejection: If a chain contains a single non-standard residue (e.g., Amino Acids, Modified Bases, or Unknown 'N'), the entire chain is discarded.

This ensures that our distance calculations are based solely on pure, unambiguous RNA interactions.

## Training the scoring function

The training script (rna_training.py) processes the cleaned dataset to generate the "scoring profiles". The goal is to learn the statistical rules of native RNA folding by following three key steps:

- Distance Categorization (The Rules) We extract geometric data to define what a "contact" looks like:

    - Atoms: We measure distances between C3' atoms.

    - Filtering: We only consider intrachain interactions (within the same molecule) where residues are separated by at least 3 positions (∣i−j∣>3). This focuses the model on folding rather than local bonds.

    - Grouping: Pairs are grouped into 10 categories (AA, AU, AC, AG, UU, UC, UG, CC, CG, GG) plus a reference category.

- Frequency Calculation (The Statistics) We calculate how often interactions occur compared to random chance:

    - Observed Frequency (fijObs​): How often do we see a specific pair (e.g., A-U) at a certain distance? 

    - Reference Frequency (fXXRef​): How often do we see any pair at that distance? This "XX" reference represents the average shape of the RNA backbone, ignoring specific base types.

- Score Computation (The Energy) We convert these probabilities into a pseudo-energy score using the Inverse Boltzmann principle. The formula is:
u(r)=−log(fXXRef​(r)fijObs​(r)​)

    - Negative Score: The pair is observed more often than the reference. This indicates a favorable, stable interaction (Low Energy).

    - Positive Score: The pair is observed less often than the reference. This indicates an unfavorable interaction (High Energy).

The script outputs 10 scoring files, which serve as the "knowledge base" for evaluating new RNA structures.

# Usage

To train the objective function, run the rna_training.py script from the terminal. You must specify the input folder and the file format.

Option 1: Training with PDB files

If your data is in rna_data/pdb:
Bash

py rna_training.py --data rna_data/pdb --format pdb

Option 2: Training with mmCIF files

If your data is in rna_data/cif:
Bash

py rna_training.py --data rna_data/cif --format cif

Arguments

    -d or --data: Path to the directory containing your structure files.

    -f or --format: The file extension to look for (pdb or cif). Default is pdb.

How to run it right now (in your terminal)

Since you are on Windows and had that environment issue earlier, always use py.

1. To train on PDBs:
PowerShell

py rna_training.py -d rna_data/pdb -f pdb

2. To train on CIFs:
PowerShell

py rna_training.py -d rna_data/cif -f cif


