# rna_downloader.py
# A script to use the PDB REST API to download RNA structure files
# PDB and mmCIF formats
# python rna_downloader.py

import requests
import os
import time


class PDBRNADownloader:
    def __init__(self, output_dir="data", file_format="cif", organism="Homo sapiens"):
        self.output_dir = output_dir
        self.format = file_format
        self.organism = organism

        os.makedirs(self.output_dir, exist_ok=True)
        
        if self.format not in ['pdb', 'cif']:
            raise ValueError("format must be 'pdb' or 'cif'")
    
    def search_rna(self, num_rna=20):
        """
        Query a number of RNA structures using PDB REST API
        
        Args:
            num_rna: Number of RNA structures to retrieve
            
        Returns:
            List of PDB IDs
        """
        # PDB REST API query
        query = {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name",
                            "operator": "exact_match",
                            "value": self.organism
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "entity_poly.rcsb_entity_polymer_type",
                            "operator": "exact_match",
                            "value": "RNA"
                        }
                    }
                ]
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {
                    "start": 0,
                    "rows": num_rna
                }
            }
        }
        
        url = "https://search.rcsb.org/rcsbsearch/v2/query"
        
        try:
            response = requests.post(url, json=query, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            pdb_ids = [item['identifier'] for item in data.get('result_set', [])]
            print(f"Found {len(pdb_ids)} RNA structures")
            return pdb_ids
            
        except requests.exceptions.RequestException as e:
            print(e)
            return []
    
    def download_single(self, pdb_id):
        """
        Download one structure file
        
        Args:
            pdb_id: PDB identifier (e.g., '1EHZ')
            
        Returns:
            Path to downloaded file or None if failed
        """
        pdb_id = pdb_id.upper()
        
        url = f"https://files.rcsb.org/download/{pdb_id}.{self.format}"
        filename = f"{pdb_id}.{self.format}"
    
        filepath = os.path.join(self.output_dir, filename)
        
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            with open(filepath, 'w') as f:
                f.write(response.text)
            
            print(f" {pdb_id}: Downloaded successfully")
            return filepath
            
        except requests.exceptions.RequestException as e:
            print(f" {pdb_id}: Failed to download - {e}")
            return None
    
    def download_multiple(self, pdb_ids, delay=0.5):
        """
        Download multiple structures
        
        Args:
            pdb_ids: List of PDB IDs
            
        Returns:
            List of successfully downloaded file paths
        """
        downloaded = []
        
        print(f"\nDownloading {len(pdb_ids)} structures as {self.format.upper()} format...")
        print(f"Output directory: {self.output_dir}\n")
        
        for i, pdb_id in enumerate(pdb_ids, 1):
            print(f"[{i}/{len(pdb_ids)}]", end=" ")
            filepath = self.download_single(pdb_id)
            
            if filepath:
                downloaded.append(filepath)
            
            if i < len(pdb_ids):
                time.sleep(delay)
        
        print(f"\nDownload completed for {len(downloaded)}/{len(pdb_ids)}")
        return downloaded
    
    def download_from_file(self, id_file, delay=0.5):
        """
        Download structures from a list of PDB IDs in a text file
        
        Args:
            id_file: Path to file with one PDB ID per line
            delay: Delay between downloads in seconds
        """
        with open(id_file, 'r') as f:
            pdb_ids = [line.strip() for line in f if line.strip()]
        return self.download_multiple(pdb_ids, delay)
  

if __name__ == "__main__":
    downloader = PDBRNADownloader(
        output_dir="data/pdb",
        file_format="pdb" 
    )
    
    # Random structures
    # pdb_ids = downloader.search_rna(num_rna=50)
    # if pdb_ids:
    #     downloaded_files = downloader.download_multiple(pdb_ids)
    #     print(f"\nFiles saved to: {os.path.abspath(downloader.output_dir)}")
    
    # Read the structures to download from a text file
    # !!! RUN THIS FIRST:
    # ls -1 rna_data/pdb | sed -e 's/\..*$//' > data20.txt
    downloader_cif = PDBRNADownloader(
        output_dir="data/cif",
        file_format="cif" 
    )
    downloader_cif.download_from_file('data50.txt')
