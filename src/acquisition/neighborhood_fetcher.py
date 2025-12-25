import os
import time
from dotenv import load_dotenv
from Bio import Entrez, SeqIO
from tqdm import tqdm

load_dotenv()

# Always required by NCBI
Entrez.email = os.getenv("NCBI_EMAIL") 

# This will be None if you didn't put it in .env, 
# and Biopython will automatically skip sending it.
Entrez.api_key = os.getenv("NCBI_API_KEY") 

if not Entrez.api_key:
    print("Running in Anonymous Mode (3 requests/sec limit)")
else:
    print("API Key detected (10 requests/sec limit)")

RAW_DIR = "data/raw"
INTERIM_DIR = "data/interim"

def search_homologs(anchor_id, limit=10):
    """Finds homologs of the anchor protein using NCBI BLAST."""
    print(f"Finding top {limit} homologs for {anchor_id}...")
    # Using E-search to find similar proteins (simplified BLAST-like approach)
    handle = Entrez.esearch(db="protein", term=f"{anchor_id}", retmax=limit)
    record = Entrez.read(handle)
    return record["IdList"]

def fetch_neighborhood(protein_id, window=10000):
    """Fetches genomic context. Implements local caching."""
    cache_path = os.path.join(RAW_DIR, f"{protein_id}_context.gbk")
    
    if os.path.exists(cache_path):
        return SeqIO.read(cache_path, "genbank")

    try:
        # Get IPG report to find the genomic coordinates
        handle = Entrez.efetch(db="protein", id=protein_id, rettype="ipg", retmode="xml")
        ipg_data = Entrez.read(handle)
        
        # Extract location (using first coordinates found)
        loc = ipg_data['IPGReport']['ProteinList'][0]['CDSList'][0]
        acc = loc.attributes['accver']
        start = int(loc.attributes['start'])
        stop = int(loc.attributes['stop'])
        
        # Define window and fetch GenBank
        s_start, s_stop = max(0, start - window), stop + window
        handle = Entrez.efetch(db="nucleotide", id=acc, seq_start=s_start, 
                               seq_stop=s_stop, rettype="gbwithparts", retmode="text")
        
        record = SeqIO.read(handle, "genbank")
        with open(cache_path, "w") as f:
            SeqIO.write(record, f, "genbank")
        
        time.sleep(0.5) # Mandatory rate limit to be a good citizen
        return record
    except Exception as e:
        print(f" Error fetching {protein_id}: {e}")
        return None

def run_pilot(anchor_id):
    homologs = search_homologs(anchor_id)
    neighborhoods = []
    
    for hit in tqdm(homologs, desc="Fetching Neighborhoods"):
        nb = fetch_neighborhood(hit)
        if nb:
            neighborhoods.append(nb)
    
    print(f"✅ Successfully cached {len(neighborhoods)} neighborhoods.")
    return neighborhoods

if __name__ == "__main__":
    # Use your Methyl Parathion Hydrolase ID as the test anchor
    test_anchor = "WP_011019623.1" 
    run_pilot(test_anchor)