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

def get_hit_details(id_list):
    """Fetches and prints titles for the IDs found by esearch."""
    if not id_list:
        print("No hits found to summarize.")
        return
    
    print(f"\n--- Identifying {len(id_list)} Hits ---")
    # esummary gives us the 'Title' and 'Organism' without downloading the whole genome
    handle = Entrez.esummary(db="protein", id=",".join(id_list))
    summaries = Entrez.read(handle)
    handle.close()
    
    for i, entry in enumerate(summaries):
        # entry['Title'] is usually: "Protein Name [Organism Name]"
        print(f"Hit {i+1}: {entry['Id']} - {entry['Title']}")


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

# Pipeline starts here
def run_pilot(anchor_id):
    homologs = search_homologs(anchor_id)

    get_hit_details(homologs)

    neighborhoods = []
    
    for hit in tqdm(homologs, desc="Fetching Neighborhoods"):
        nb = fetch_neighborhood(hit)
        if nb:
            neighborhoods.append(nb)
    
    print(f"✅ Successfully cached {len(neighborhoods)} neighborhoods.")
    return neighborhoods


def print_neighborhood_summary(neighborhood_record):
    """Prints a summary of genes found in the fetched genomic slice."""
    print(f"\n--- Neighborhood for: {neighborhood_record.annotations.get('organism', 'Unknown')} ---")
    print(f"Source Accession: {neighborhood_record.id}")
    
    # Filter for 'CDS' features (Coding Sequences)
    genes = [f for f in neighborhood_record.features if f.type == "CDS"]
    
    print(f"{'Gene/Protein ID':<20} {'Start':<10} {'End':<10} {'Product'}")
    print("-" * 70)
    
    for gene in genes:
        # Extract common qualifiers
        locus_tag = gene.qualifiers.get("locus_tag", ["N/A"])[0]
        protein_id = gene.qualifiers.get("protein_id", ["N/A"])[0]
        product = gene.qualifiers.get("product", ["N/A"])[0]
        
        start = gene.location.start
        end = gene.location.end
        
        print(f"{protein_id:<20} {start:<10} {end:<10} {product}")

# Update your __main__ block to use it:
if __name__ == "__main__":
    test_query = "Methyl parathion hydrolase"
    neighborhoods = run_pilot(test_query)
    
    for nb in neighborhoods[:2]:
        print_neighborhood_summary(nb)