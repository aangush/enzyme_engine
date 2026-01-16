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
    for feat in record.features:
        if feat.type == "CDS":
            # The coordinates in the 'record' are relative to the slice (0 to ~20000)
            # The anchor is always roughly at index 'window' (10000)
            f_start = int(feat.location.start)
            
            # Skip if this is the anchor itself
            if f_start >= (window - 100) and f_start <= (window + 100):
                continue
                
            product = feat.qualifiers.get("product", ["unknown"])[0]
            prot_id = feat.qualifiers.get("protein_id", ["no_id"])[0]
            translation = feat.qualifiers.get("translation", [""])[0]
            dist = int(abs(window - f_start))
            
            db.add_neighbor(instance_id, prot_id, product, dist, "N/A", translation)

def generate_summary_report():
    print("\n" + "="*55)
    print("📊 GNN LINKAGE SUMMARY")
    print("="*55)
    conn = sqlite3.connect("data/scout.db")
    cursor = conn.cursor()
    cursor.execute("SELECT product, COUNT(*) as freq FROM neighbors GROUP BY product ORDER BY freq DESC LIMIT 15")
    for row in cursor.fetchall():
        print(f"{row[0][:40]:<40} | {row[1]}")
    conn.close()
    print("="*55)

if __name__ == "__main__":
    run_gnn_scout("input.fasta", hit_limit=15)