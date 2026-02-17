import os
import sys
import sqlite3
import time

from discovery import DiscoveryEngine
from ncbi_client import NCBIClient
from db_manager import GNNDB
from domain_analyzer import DomainAnalyzer

def run_gnn_scout(input_fasta_path, hit_limit=200):
    discovery = DiscoveryEngine()
    client = NCBIClient()
    db = GNNDB()

    with open(input_fasta_path, "r") as f:
        fasta_str = f.read()

    # 1. Initial BLAST
    homolog_ids, identities = discovery.find_homologs(fasta_str, hit_limit=hit_limit)

    # 2. GNN SEARCHING
    print(f"Starting search for {len(homolog_ids)} homologs. (Press Ctrl+C once to skip a stuck ID)")
    for i, pid in enumerate(homolog_ids, 1):
        try:
            # Catches Ctrl+C to allow manual skipping of stalled network calls
            try:
                print(f"[{i}/{len(homolog_ids)}] Processing {pid}...")
                
                locations = client.get_locations(pid)
                if not locations: continue
            
                loc = locations[0]
                loc['anchor_id'] = pid
            
                if db.instance_exists(pid, loc['nuc_acc']):
                    print(f"Already in database.")
                    continue

                record = client.fetch_neighborhood(loc['nuc_acc'], loc['start'], loc['end'])
            
                if record:
                    instance_id = db.add_instance(pid, loc['nuc_acc'], loc['start'], loc['end'], loc['strand'])
                    extract_and_save_neighbors(record, instance_id, loc, db)
                    print(f"Neighborhood for {pid} saved.")

            except KeyboardInterrupt:
                print(f"\nManual skip triggered for {pid}. Moving to next...")
                time.sleep(1) # Brief pause to clear the signal buffer
                continue

        except Exception as e:
            print(f"Unexpected error on {pid}: {e}")
            continue

    # 3. DOMAIN ANALYSIS
    print("\nStarting Post-Processing: Domain Analysis...")
    analyzer = DomainAnalyzer(db_path="data/GNN.db")
    analyzer.analyze_hypotheticals()

    # 3.5. Display Initial BLAST Summary
    # Display BLAST Distribution
    if identities:
        display_blast_metrics(identities)

    # 4. REPORT
    generate_summary_report()
    


def display_blast_metrics(identities):
    # Calculate and display basic distribution of BLAST % Identities

    hi = max(identities)
    lo = min(identities)
    avg = sum(identities) / len(identities)

    print("\n" + "-"*60)
    print(f"Initial BLAST Info (n={len(identities)})")
    print(f"Range: {lo:.1f}% - {hi:.1f}% Identity | Average: {avg:.1f}%")

    # ASCII Distribution (Binned by 10%)
    bins = {90:0, 80:0, 70:0, 60:0, 50:0, 40:0, 30:0, 20:0}
    for ident in identities:
        for b in sorted(bins.keys(), reverse=True):
            if ident >= b:
                bins[b] += 1
                break

    print("\n% Identity Distribution:")
    for b, count in sorted(bins.items(), reverse=True):
        bar = "█" * int(count / (len(identities)/20)) if count > 0 else ""
        print(f"  {b}%+: {count:<4} {bar}")
    print("-" * 60 + "\n")

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

def suggest_module(co_results):
    # Build a map of protein neighbors from co-occurrence data
    neighbors_map = {}

    for a, b, _ in co_results:
        if a not in neighbors_map:
            neighbors_map[a] = set()
        if b not in neighbors_map:
            neighbors_map[b] = set()

        neighbors_map[a].add(b)
        neighbors_map[b].add(a)

    # Walk map to find connected groups
    visited = set()
    all_modules = []

    for protein in neighbors_map:
        if protein not in visited:
            current_module = set()
            stack = [protein]

            while stack:
                node = stack.pop()
                if node not in visited:
                    visited.add(node)
                    current_module.add(node)
                    stack.extend(neighbors_map - visited)

            all_modules.append(current_module)

    return all_modules
        
    
def generate_summary_report():
    print("\n" + "="*145)
    print(f"{'GENOMIC NEIGHBORHOOD & CO-OCCURRENCE REPORT':^145}")
    print("="*145)
    
    conn = sqlite3.connect("data/GNN.db")
    cursor = conn.cursor()
    
    # 1. Main Neighbor Report with Representative Accession
    # We use MAX(protein_id) to just grab one valid accession as an example
    query = """
    SELECT 
        product, 
        MAX(protein_id) as example_acc,
        direction, 
        COUNT(*) as freq, 
        CAST(AVG(distance_bp) AS INT) as avg_dist,
        CAST(MAX(distance_bp) - MIN(distance_bp) AS INT) as spread,
        CAST(AVG(LENGTH(sequence)) AS INT) as avg_len
    FROM neighbors 
    GROUP BY product, direction
    HAVING freq >= 2
    ORDER BY freq DESC 
    LIMIT 40
    """
    cursor.execute(query)
    results = cursor.fetchall()
    
    # Header with expanded padding for Name + Accession
    print(f"{'Neighbor Product (Example Accession)':<65} | {'Strand':<8} | {'Freq':<5} | {'Avg Dist':<12} | {'Spread':<8} | {'Len'}")
    print("-" * 145)
    
    for row in results:
        # Combine Product Name and Accession for the display string
        display_name = f"{row[0]} ({row[1]})"
        # Truncate if too long for the 65-char column
        display_name = (display_name[:62] + '..') if len(display_name) > 65 else display_name
        
        print(f"{display_name:<65} | {row[2]:<8} | {row[3]:<5} | {row[4]:>9} bp | {row[5]:>6} bp | {row[6]}aa")

    # 2. (Co-occurrence) Section
    print("\n" + "-"*60)
    print("TOP CO-OCCURRING NEIGHBOR PAIRS (Pathway Modules)")
    print("-"*60)
    
    co_query = """
    SELECT 
        n1.product AS Neighbor_A, 
        n2.product AS Neighbor_B, 
        COUNT(DISTINCT n1.instance_id) AS CoOccurrence_Freq
    FROM neighbors n1
    JOIN neighbors n2 ON n1.instance_id = n2.instance_id
    WHERE n1.product < n2.product 
    GROUP BY Neighbor_A, Neighbor_B
    HAVING CoOccurrence_Freq >= 3
    ORDER BY CoOccurrence_Freq DESC
    LIMIT 15
    """
    cursor.execute(co_query)
    co_results = cursor.fetchall()
    
    if co_results:
        print(f"{'Module Component A':<40} | {'Module Component B':<40} | {'Freq'}")
        for crow in co_results:
            print(f"{crow[0][:39]:<40} | {crow[1][:39]:<40} | {crow[2]}")
    else:
        print("No significant co-occurrence modules found yet.")
        
    conn.close()
    print("="*145)

    modules = suggest_module(co_results)

    if modules:
        print("\n--- Putative Functional Modules (Operons) ---")
        for i, mod in enumerate(modules, 1):
        # Joining the set into a readable string
            print(f"Module {i}: {', '.join(mod)}")

if __name__ == "__main__":
    run_gnn_scout("input.fasta", hit_limit=200)
    
