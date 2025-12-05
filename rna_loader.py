import os
import sys
from Bio.PDB import PDBParser, MMCIFParser

def load_rna_structure(file_path):
    """
    Loads a PDB or mmCIF file and returns the first model.
    
    Args:
        file_path (str): Path to the .pdb or .cif file.
        
    Returns:
        Bio.PDB.Model.Model: The first model of the structure, or None if failed.
    """
    # 1. Detect file format
    _, ext = os.path.splitext(file_path)
    ext = ext.lower()

    # 2. Initialize the correct parser
    if ext == '.pdb':
        parser = PDBParser(QUIET=True)
    elif ext in ['.cif', '.mcif']:
        parser = MMCIFParser(QUIET=True)
    else:
        print(f"Error: Unsupported format '{ext}' for file {file_path}")
        return None

    try:
        # 3. Parse and return the first model (Model 0)
        structure_id = os.path.basename(file_path)
        structure = parser.get_structure(structure_id, file_path)
        return structure[0]
        
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None