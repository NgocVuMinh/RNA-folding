# RNA statistical potentials: knowledge-based scoring function

**Course: Bioinformatics of RNA and non-coding world**

**M2 GENIOMHE - Univ Evry Paris-Saclay**

**Year: 2025-2026**

**Authors:**  
- Minh Ngoc VU  
- Erine BENOIST

---

## Overview

For a given ribonucleotide chain, the RNA folding problem consists of finding the native fold among an astronomically large number of possible conformations. The native fold is generally considered the one with the lowest Gibbs free energy.

The goal of this project is to develop an objective function (scoring function) to estimate this energy. By evaluating the “pseudo-energy” of a conformation, we can determine how close a predicted structure is to the optimal/native fold.

Our scoring function is a statistical potential derived from experimentally determined RNA structures. It relies on the inverse Boltzmann principle, calculating pseudo-energies based on the frequency of atomic interactions observed in nature.

---

## Installation & Requirements

The project was developed using python 3.12.

### 1. Dependencies

Install the required packages using pip:

```bash
pip install biopython numpy scipy matplotlib
```

### 2. Project structure

- **main.py** — Central entry point, handles CLI arguments and orchestrates the training process  
- **rna_loader.py** — Parses PDB and CIF files using BioPython  
- **rna_distance.py** — Computes Euclidean distances and filters invalid/non-RNA chains  
- **rna_training.py** — Core statistical potential logic (Histograms + KDE)  
- **rna_plot.py** — Plot profiles
- **rna_downloader.py** — Downloads structures via the RCSB PDB API  
- **utils.py** — Helper functions (file saving, formatting)

---

## Usage

### 1. Data acquisition

RNA structures can be downloaded using:

```bash
python rna_downloader.py
```
This will create an `rna_data/` folder containing PDB and mmCIF files.

For training, we obtained 114 RNA structures from *Homo sapiens* via the [RCSB PDB REST API](https://www.rcsb.org/docs/programmatic-access/web-apis-overview). Each structure was downloaded in both PDB and CIF format.

---

### 2. Training the objective function

#### Basic usage

Using histograms:

```bash
python train.py --data data/pdb --format pdb --out_dir potentials
```

#### Advanced usage

Kernel Density Estimation is still being developed. Users can specify atom of choice and mmCIF format. The project currently supports single-atom distances.

```bash
python train.py --data data/cif --format cif --mode kernel --atom "C3'" --out_dir out_kde
```

| Argument | Description | Default |
|---------|-------------|---------|
| `-d, --data` | Path to structure folder | **Required** |
| `-f, --format` | File format (`pdb` or `cif`) | pdb |
| `-m, --mode` | Calculation method (`histogram` or `kernel`) | histogram |
| `-a, --atom` | Atom used for distance calculation | C3' |
| `--max_dist` | Maximum distance threshold (Å) | 20.0 |
| `-o, --out_dir` | Output directory | potentials |

---

### 3. Plot interaction profiles

```bash
python plot.py --input out_hist --output plots_hist
```

Parameters:  
- `--input` : Folder containing .txt scoring files
- `--output`: Folder where .png plots will be saved  

---

### 4. Scoring on a given structure 

Scoring based on histograms are supported. KDE is still being developed.

```bash
python scoring.py --structure_dir data/pdb/1A3M.pdb --profile_dir out_hist
```

Parameters:  
- `--structure_dir`: RNA structure file in PDB or mmCIF
- `--profile_dir`: path to where the profiles are stored 


---

## Methodology

### 1. Data loading & cleaning

- Uses **Biopython** to parse structures  
- Handles **both PDB and CIF formats**  
- Selects **Model 0** for NMR files  
- Removes chains containing:
  - amino acids  
  - DNA  
  - modified bases  
  - ambiguous nucleotide “N”  

Only **pure RNA chains** are kept for training.

---

### 2. Distance calculation

We define an “interaction” as:

- Measured between **C3' atoms**  
- Only **intrachain** interactions  
- Residue separation ≥ 3 positions:  

\[
|i - j| > 3
\]

This avoids trivial backbone-local interactions.

---

### 3. The objective function

The statistical potential is computed using the inverse Boltzmann principle:

\[
u(r) = -\log\left(\frac{f_{ij}^{Obs}(r)}{f_{XX}^{Ref}(r)}\right)
\]

Where:

- \( f_{ij}^{Obs}(r) \): Observed probability of interacting pair \(i-j\) at distance \(r\)  
- \( f_{XX}^{Ref}(r) \): Reference state probability (any pair at distance \(r\))  

**Interpretation of the score:**

- **Negative** → interaction more frequent than random → **favorable**  
- **Positive** → interaction less frequent than random → **unfavorable**

