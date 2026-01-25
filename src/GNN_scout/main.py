import os
import sys
import sqlite3
import time

# Ensure local imports work
sys.path.append(os.path.dirname(__file__))

from discovery import DiscoveryEngine
from ncbi_client import NCBIClient
from db_manager import ScoutDB
from domain_analyzer import DomainAnalyzer

def run_gnn_scout(input_fasta_path, hit_limit=100):
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

    # 2. SEARCHING
    print(f"🚀 Starting scout for {len(homolog_ids)} homologs. (Press Ctrl+C once to skip a stuck ID)")
    for i, pid in enumerate(homolog_ids, 1):
        try:
            # INNER SHIELD: Catches Ctrl+C to allow manual skipping of stalled network calls
            try:
                print(f"[{i}/{len(homolog_ids)}] Processing {pid}...")
                
                locations = client.get_locations(pid)
                if not locations: continue
            
                loc = locations[0]
                loc['anchor_id'] = pid 
            
                if db.instance_exists(pid, loc['nuc_acc']):
                    print(f"⏩ Already in database.")
                    continue

                record = client.fetch_neighborhood(loc['nuc_acc'], loc['start'], loc['end'])
            
                if record:
                    instance_id = db.add_instance(pid, loc['nuc_acc'], loc['start'], loc['end'], loc['strand'])
                    extract_and_save_neighbors(record, instance_id, loc, db)
                    print(f"✅ Neighborhood for {pid} saved.")

            except KeyboardInterrupt:
                print(f"\n🛑 Manual skip triggered for {pid}. Moving to next...")
                time.sleep(1) # Brief pause to clear the signal buffer
                continue

        except Exception as e:
            print(f"❌ Unexpected error on {pid}: {e}")
            continue

    # 3. DOMAIN ANALYSIS
    print("\n🕵️ Starting Post-Processing: Domain Analysis...")
    analyzer = DomainAnalyzer(db_path="data/scout.db")
    analyzer.analyze_hypotheticals()

    # 4. REPORT
    generate_summary_report()

def extract_and_save_neighbors(record, instance_id, loc, db, window=10000):
    anchor_id = loc.get('anchor_id')
    anchor_strand = loc['strand']
    
    for feat in record.features:
        if feat.type == "CDS":
            prot_id = feat.qualifiers.get("protein_id", ["no_id"])[0]
            if prot_id == anchor_id: continue
                
            f_start = int(feat.location.start)
            f_end = int(feat.location.end)
            
            # Simple spatial exclusion to prevent self-overlap issues
            if f_start <= window <= f_end: continue

            product = feat.qualifiers.get("product", ["unknown"])[0]
            translation = feat.qualifiers.get("translation", [""])[0]
            
            # Strand Calculation
            neighbor_strand = feat.location.strand
            direction = "Same" if neighbor_strand == anchor_strand else "Oppose"
            
            # Signed Distance Calculation (Center of neighbor - Anchor position)
            feat_center = (f_start + f_end) / 2
            dist = int(feat_center - window)
            
            db.add_neighbor(instance_id, prot_id, product, dist, direction, translation)

def generate_summary_report():
    print("\n" + "="*115)
    print("📊 GENOMIC NEIGHBORHOOD REPORT: SYNTENY & DIVERSITY MAPPING")
    print("="*115)
    
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
    
    print(f"{'Neighbor Product':<50} | {'Strand':<8} | {'Freq':<5} | {'Avg Dist':<12} | {'Spread':<8} | {'Avg Len'}")
    print("-" * 115)
    for row in results:
        p_name = (row[0][:47] + '..') if len(row[0]) > 47 else row[0]
        print(f"{p_name:<50} | {row[1]:<8} | {row[2]:<5} | {row[3]:>9} bp | {row[4]:>6} bp | {row[5]}aa")
    
    conn.close()
    print("="*115)

if __name__ == "__main__":
    run_gnn_scout("input.fasta", hit_limit=100)