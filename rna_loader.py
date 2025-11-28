# rna_loader.py
# I used 1EHZ as reference for RNA structure loading
# Dowload the files using:
# curl https://files.rcsb.org/download/1EHZ.pdb -o 1EHZ.pdb
# curl https://files.rcsb.org/download/1EHZ.cif -o 1EHZ.cif

import os
from Bio.PDB import PDBParser, MMCIFParser

def load_rna_structure(file_path):
    """
    Loads a PDB or mmCIF file and returns the first model.
    """
    # Get the file extension
    _, ext = os.path.splitext(file_path)
    ext = ext.lower()

    # Select the correct parser
    if ext == '.pdb':
        parser = PDBParser(QUIET=True)
    elif ext in ['.cif', '.mcif']:
        parser = MMCIFParser(QUIET=True)
    else:
        raise ValueError(f"Format {ext} not supported. Use .pdb or .cif")

    try:
        # Parse the structure
        structure_id = os.path.basename(file_path)
        structure = parser.get_structure(structure_id, file_path)
        
        # Return the first model (standard for RNA/NMR)
        return structure[0]
        
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None