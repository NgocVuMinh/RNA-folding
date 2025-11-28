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

    1. Standard Bases Only: It scans every chain to ensure it is composed exclusively of standard RNA bases (Adenine, Uracil, Cytosine, Guanine).

    2. Whole-Chain Rejection: If a chain contains a single non-standard residue (e.g., Amino Acids, Modified Bases, or Unknown 'N'), the entire chain is discarded.

This ensures that our distance calculations are based solely on pure, unambiguous RNA interactions.

## Training the scroring function

The training script (rna_training.py) processes the cleaned dataset to generate the "scoring profiles" (objective function). The process involves three key steps:

    1. Distance Categorization:

    - It extracts C3' atoms and calculates Euclidean distances between residue pairs.
    - Only intrachain distances are considered for training, with a sequence separation of at least 3 residues (∣i−j∣>3).
    - Pairs are grouped into 10 categories (AA, AU, AC, AG, UU, UC, UG, CC, CG, GG), where order is ignored (e.g., AU and UA are counted together).

    2. Frequency Calculation:

    - Observed Frequency (fijObs​): The probability of observing a specific pair (e.g., A-U) at a specific distance bin r.
    - Reference Frequency (fXXRef​): A "virtual" reference state calculated by counting all valid pairs regardless of residue type.

     3. Score Computation: The final pseudo-energy score for each distance bin is computed using the negative log-likelihood ratio:

            u(r)=−log(fXXRef​(r)fijObs​(r)​)

The script outputs 10 scoring files, representing the energy profiles for each base pair type.



