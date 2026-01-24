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
        loc['anchor_id'] = pid

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
    # Anchor metadata for comparison
    anchor_id = loc.get('anchor_id')
    anchor_strand = loc['strand']
    
    for feat in record.features:
        if feat.type == "CDS":
            prot_id = feat.qualifiers.get("protein_id", ["no_id"])[0]
            
            # 1. Identity Check: Skip if this neighbor is actually the anchor protein
            if prot_id == anchor_id:
                continue
                
            # 2. Coordinate Check: Skip if it overlaps the center of our window
            f_start = int(feat.location.start)
            f_end = int(feat.location.end)
            if f_start <= window <= f_end:
                continue

            product = feat.qualifiers.get("product", ["unknown"])[0]
            translation = feat.qualifiers.get("translation", [""])[0]
            
            # 3. Dynamic Strand Calculation (The Fix)
            # Biopython uses 1 for forward, -1 for reverse
            neighbor_strand = feat.location.strand
            direction = "Same" if neighbor_strand == anchor_strand else "Oppose"
            
            # 4. Distance Calculation
            feat_center = (f_start + f_end) / 2
            dist = int(abs(window - feat_center))
            
            # Save with the real direction instead of "N/A"
            db.add_neighbor(instance_id, prot_id, product, dist, direction, translation)

def generate_summary_report():
    print("\n" + "="*75)
    print("📊 GNN LINKAGE SUMMARY: RELATIVE STRAND ANALYSIS")
    print("="*75)
    
    db_path = "data/scout.db"
    if not os.path.exists(db_path):
        print("❌ Database not found.")
        return

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Query: Group by product and direction to see conservation patterns
    query = """
    SELECT product, direction, COUNT(*) as freq 
    FROM neighbors 
    GROUP BY product, direction 
    ORDER BY freq DESC 
    LIMIT 25
    """
    cursor.execute(query)
    results = cursor.fetchall()
    
    print(f"{'Neighbor Product':<45} | {'Strand':<10} | {'Freq'}")
    print("-" * 75)
    for row in results:
        # Format the display name
        display_name = (row[0][:42] + '..') if len(row[0]) > 42 else row[0]
        # row[1] is now the actual 'Same' or 'Oppose' string
        print(f"{display_name:<45} | {row[1]:<10} | {row[2]}")
    
    conn.close()
    print("="*75)

if __name__ == "__main__":
    run_gnn_scout("input.fasta", hit_limit=20)