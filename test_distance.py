from rna_loader import load_rna_structure
from rna_distance import get_all_distances
import os

if __name__ == "__main__":
    # Test with your file
    test_file = "1EHZ.pdb" 
    
    if os.path.exists(test_file):
        print(f"Processing {test_file}...")
        model = load_rna_structure(test_file)
        
        # Calculate!
        data = get_all_distances(model)
        
        # Show results
        print(f"Found {len(data)} valid interactions (< 20 A).")
        
        # Count types
        intra = len([x for x in data if x['Type'] == 'Intrachain'])
        inter = len([x for x in data if x['Type'] == 'Interchain'])
        
        print(f"Intrachain: {intra}")
        print(f"Interchain: {inter}")
        
        # Show a sample
        if len(data) > 0:
            print("\nSample Interaction:")
            print(data[0])