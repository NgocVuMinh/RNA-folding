from rna_loader import load_rna_structure

# Test PDB
pdb_model = load_rna_structure("1X8W.pdb")
print(pdb_model)
print(f"PDB Loaded: {pdb_model.id} | Chains: {[c.id for c in pdb_model]}")

# Test mmCIF
cif_model = load_rna_structure("1X8W.cif")
print(f"CIF Loaded: {cif_model.id} | Chains: {[c.id for c in cif_model]}")