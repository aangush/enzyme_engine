import os
import sys
import sqlite3

# Ensure local imports work
sys.path.append(os.path.dirname(__file__))

from discovery import DiscoveryEngine
from ncbi_client import NCBIClient
from db_manager import ScoutDB
from domain_analyzer import DomainAnalyzer

def run_gnn_scout(input_fasta_path, hit_limit=50):
    discovery = DiscoveryEngine()
    client = NCBIClient()
    db = ScoutDB()
    
    if not os.path.exists(input_fasta_path):
        print(f"Input file {input_fasta_path} not found.")
        return

    with open(input_fasta_path, "r") as f:
        fasta_str = f.read()

    # 1. DISCOVERY - Sampling the first 50 hits
    homolog_ids = discovery.find_homologs(fasta_str, hit_limit=hit_limit)

    # 2. SCOUTING
    for pid in homolog_ids:
        # Avoid redundant downloads
        locations = client.get_locations(pid)
        if not locations: continue
        
        loc = locations[0]
        loc['anchor_id'] = pid 
        
        if db.instance_exists(pid, loc['nuc_acc']):
            print(f"⏩ Skipping {pid}, already in database.")
            continue

        print(f"Scouting neighborhood for homolog: {pid}")
        record = client.fetch_neighborhood(loc['nuc_acc'], loc['start'], loc['end'])
        
        if record:
            instance_id = db.add_instance(pid, loc['nuc_acc'], loc['start'], loc['end'], loc['strand'])
            extract_and_save_neighbors(record, instance_id, loc, db)

    # 2.5 DOMAIN ANALYSIS
    analyzer = DomainAnalyzer(db_path="data/scout.db")
    analyzer.analyze_hypotheticals()

    # 3. REPORT - The Linkage Summary
    generate_summary_report()

def extract_and_save_neighbors(record, instance_id, loc, db, window=10000):
    anchor_id = loc.get('anchor_id')
    anchor_strand = loc['strand']
    
    for feat in record.features:
        if feat.type == "CDS":
            prot_id = feat.qualifiers.get("protein_id", ["no_id"])[0]
            
            # Identity Check
            if prot_id == anchor_id:
                continue
                
            # Spatial Check
            f_start = int(feat.location.start)
            f_end = int(feat.location.end)
            if f_start <= window <= f_end:
                continue

            product = feat.qualifiers.get("product", ["unknown"])[0]
            translation = feat.qualifiers.get("translation", [""])[0]
            
            # Strand Calculation
            neighbor_strand = feat.location.strand
            direction = "Same" if neighbor_strand == anchor_strand else "Oppose"
            
            # Distance
            feat_center = (f_start + f_end) / 2
            dist = int(abs(window - feat_center))
            
            db.add_neighbor(instance_id, prot_id, product, dist, direction, translation)

def generate_summary_report():
    print("\n" + "="*95)
    print("GNN ANALYSIS: TOP CONSERVED GENOMIC NEIGHBORS")
    print("="*95)
    
    conn = sqlite3.connect("data/scout.db")
    cursor = conn.cursor()
    
    # We only show neighbors appearing multiple times (High-frequency linkage)
    query = """
    SELECT product, direction, COUNT(*) as freq, CAST(AVG(distance_bp) AS INT) as avg_dist
    FROM neighbors 
    GROUP BY product, direction 
    HAVING freq >= 1
    ORDER BY freq DESC 
    LIMIT 30
    """
    cursor.execute(query)
    results = cursor.fetchall()
    
    print(f"{'Neighbor Product (Identified Domains)':<50} | {'Strand':<10} | {'Freq'} | {'Avg Dist'}")
    print("-" * 80)
    for row in results:
        p_name = (row[0][:47] + '..') if len(row[0]) > 47 else row[0]
        print(f"{p_name:<50} | {row[1]:<10} | {row[2]:<5} | {row[3]} bp")
    
    conn.close()
    print("="*95)

if __name__ == "__main__":
    # Target 50 homologs to map the diversity space
    run_gnn_scout("input.fasta", hit_limit=50)