import os
import sys
import sqlite3

# Ensure local imports work regardless of run location
sys.path.append(os.path.dirname(__file__))

from discovery import DiscoveryEngine
from ncbi_client import NCBIClient
from db_manager import ScoutDB

def run_gnn_scout(input_fasta_path, hit_limit=15):
    discovery = DiscoveryEngine()
    client = NCBIClient()
    db = ScoutDB()
    
    if not os.path.exists(input_fasta_path):
        print(f"❌ Input file {input_fasta_path} not found.")
        return

    with open(input_fasta_path, "r") as f:
        fasta_str = f.read()

    # 1. DISCOVERY
    homolog_ids = discovery.find_homologs(fasta_str, hit_limit=hit_limit)

    # 2. SCOUTING
    for pid in homolog_ids:
        # Avoid duplicate work
        locations = client.get_locations(pid)
        if not locations: continue
        
        loc = locations[0]
        if db.instance_exists(pid, loc['nuc_acc']):
            print(f"⏩ Skipping {pid}, already in database.")
            continue

        record = client.fetch_neighborhood(loc['nuc_acc'], loc['start'], loc['end'])
        if record:
            instance_id = db.add_instance(pid, loc['nuc_acc'], loc['start'], loc['end'], loc['strand'])
            extract_and_save_neighbors(record, instance_id, loc, db)
            print(f"✅ Neighborhood for {pid} saved.")

    # 3. REPORT
    generate_summary_report()

def extract_and_save_neighbors(record, instance_id, loc, db, window=10000):
    anchor_id = loc.get('anchor_id') # Ensure this is passed or available
    anchor_strand = loc['strand']
    
    for feat in record.features:
        if feat.type == "CDS":
            prot_id = feat.qualifiers.get("protein_id", ["no_id"])[0]
            
            # 1. STRICT EXCLUSION: If the neighbor's ID is our anchor, skip it.
            if prot_id == anchor_id:
                continue
                
            # 2. COORDINATE EXCLUSION (Safety net):
            # The anchor is centered at 'window'. If the feature overlaps the center, skip.
            f_start = int(feat.location.start)
            f_end = int(feat.location.end)
            if f_start <= window <= f_end:
                continue

            product = feat.qualifiers.get("product", ["unknown"])[0]
            translation = feat.qualifiers.get("translation", [""])[0]
            
            # 3. STRAND CALCULATION:
            neighbor_strand = feat.location.strand
            direction = "Same" if neighbor_strand == anchor_strand else "Oppose"
            
            # Distance from center
            feat_center = (f_start + f_end) / 2
            dist = int(abs(window - feat_center))
            
            db.add_neighbor(instance_id, prot_id, product, dist, direction, translation)

def generate_summary_report():
    print("\n" + "="*70)
    print("📊 GNN LINKAGE SUMMARY")
    print("="*70)
    conn = sqlite3.connect("data/scout.db")
    cursor = conn.cursor()
    
    # Group by product AND direction to see if certain genes 
    # are always on the same or opposite strand as your anchor.
    query = """
    SELECT product, direction, COUNT(*) as freq 
    FROM neighbors 
    GROUP BY product, direction 
    ORDER BY freq DESC 
    LIMIT 20
    """
    cursor.execute(query)
    results = cursor.fetchall()
    
    print(f"{'Neighbor Product':<40} | {'Strand':<8} | {'Freq'}")
    print("-" * 70)
    for row in results:
        product_name = (row[0][:37] + '..') if len(row[0]) > 37 else row[0]
        print(f"{product_name:<40} | {row[1]:<8} | {row[2]}")
    conn.close()
    print("="*70)

if __name__ == "__main__":
    run_gnn_scout("input.fasta", hit_limit=15)