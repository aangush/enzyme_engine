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
            
            # Signed Distance
            feat_center = (f_start + f_end) / 2
            dist = int(feat_center - window)
            
            db.add_neighbor(instance_id, prot_id, product, dist, direction, translation)

def generate_summary_report():
    print("\n" + "="*110)
    print("GENOMIC NEIGHBORHOOD REPORT FOR INPUT FASTA: NO PFAMS GROUPED BY DISTANCE")
    print("="*110)
    
    conn = sqlite3.connect("data/scout.db")
    cursor = conn.cursor()
    
    query = """
    SELECT 
        product, 
        direction, 
        COUNT(*) as freq, 
        CAST(AVG(distance_bp) AS INT) as avg_dist,
        CAST(MAX(distance_bp) - MIN(distance_bp) AS INT) as spread,
        CAST(AVG(LENGTH(sequence)) AS INT) as avg_len
    FROM neighbors 
    GROUP BY product, direction
    HAVING freq >= 2
    ORDER BY freq DESC 
    LIMIT 100
    """
    cursor.execute(query)
    results = cursor.fetchall()
    
    print(f"{'Neighbor Product':<45} | {'Strand':<8} | {'Freq':<5} | {'Avg Dist':<10} | {'Spread':<8} | {'Avg Len'}")
    print("-" * 105)
    for row in results:
        p_name = (row[0][:42] + '..') if len(row[0]) > 42 else row[0]
        # Spread tells us if the genes are 'locked' in the same spot across different genomes
        print(f"{p_name:<50} | {row[1]:<8} | {row[2]:<5} | {row[3]:>9} bp | {row[4]:>6} bp | {int(row[5])}aa")
    
    conn.close()
    print("="*110)

if __name__ == "__main__":
    # Target 100 homologs to map the diversity space
    run_gnn_scout("input.fasta", hit_limit=100)