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

## Cleaning the chains

To generate a high-quality statistical potential, we implemented a strict filtering mechanism in `rna_distance.py`. Experimental PDB files often contain "noise" such as protein chains, DNA-RNA hybrids, or ambiguous IUPAC codes.

Our `is_strict_pure_rna` function applies the following logic:

- Standard bases only: It scans every chain to ensure it is composed exclusively of standard RNA bases (Adenine, Uracil, Cytosine, Guanine).

- Whole-chain rejection: If a chain contains a single non-standard residue (e.g., amino acids, modified bases, or unknown 'N'), the entire chain is discarded.

This ensures that our distance calculations are based solely on pure, unambiguous RNA interactions.

## Training the scoring function

The training script (`rna_training.py`) processes the cleaned dataset to generate the "scoring profiles". The goal is to learn the statistical rules of native RNA folding by following three key steps:

- Distance categorization: We extract geometric data to define what a "contact" looks like:

    - Atoms: We measure distances between a specific atom (default is C3').

    - Filtering: We only consider intrachain interactions (within the same molecule) where residues are separated by at least 3 positions (∣i−j∣>3). This focuses the model on folding rather than local bonds.

    - Grouping: Pairs are grouped into 10 categories (AA, AU, AC, AG, UU, UC, UG, CC, CG, GG) plus a reference category.

- Frequency calculation: We calculate how often interactions occur compared to random chance:

    - Observed frequency (fijObs​): How often we see a specific pair (e.g., A-U) at a certain distance.

    - Reference frequency (fXXRef​): How often we see any pair at that distance. This "XX" reference represents the average shape of the RNA backbone, ignoring specific base types.

- Score computation: We convert these probabilities into a pseudo-energy score using the Inverse Boltzmann principle. The formula is:

u(r)=−log(fXXRef​(r)fijObs​(r)​)


    - Negative score: The pair is observed more often than the reference. This indicates a favorable, stable interaction (low energy).

    - Positive score: The pair is observed less often than the reference. This indicates an unfavorable interaction (high energy).

The script outputs 10 scoring files, which serve as the "knowledge base" for evaluating new RNA structures.

## Usage

To train the objective function, run the rna_training.py script from the terminal. We specify the path to the input data and the file format (default is PDB).

Option 1: Training with PDB files

```bash
py rna_training.py --data rna_data/pdb --format pdb
```
Option 2: Training with mmCIF files

```bash
py rna_training.py --data rna_data/cif --format cif
```



