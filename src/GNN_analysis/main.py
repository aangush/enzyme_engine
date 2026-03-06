import os
import sys
import sqlite3
import asyncio
import aiohttp
from collections import defaultdict, Counter
from itertools import combinations

# Third-party imports
try:
    from tqdm.asyncio import tqdm
except ImportError:
    print("Please install tqdm: pip install tqdm")
    sys.exit(1)

# Local module imports
from discovery import DiscoveryEngine
from ena_client import UniProtENAClient
from db_manager import GNNDB

async def process_homolog(session, client, pid, window=10000):
    """
    Async worker to fetch metadata, ENA context, and return structured neighbor data.
    """
    metadata, locations = await client.fetch_uniprot_data(session, pid)
    
    if not locations or not metadata:
        return pid, None, []

    # Take the primary associated EMBL contig
    loc = locations[0]
    embl_acc = loc['embl_acc']

    record = await client.fetch_ena_neighborhood(session, embl_acc)
    if not record:
        return pid, None, []

    # ENA coordinates map differently than NCBI. Scan the record to find the anchor
    # first to establish its position and strand.
    anchor_feature = None
    for feat in record.features:
        if feat.type == "CDS":
            uniprot_refs = [ref for ref in feat.qualifiers.get("db_xref", []) if "UniProtKB" in ref]
            if any(pid in ref for ref in uniprot_refs):
                anchor_feature = feat
                break

    if not anchor_feature:
        return pid, None, []

    a_start = int(anchor_feature.location.start)
    a_end = int(anchor_feature.location.end)
    a_strand = anchor_feature.location.strand

    instance_data = {
        "embl_acc": embl_acc,
        "start": a_start,
        "end": a_end,
        "strand": a_strand
    }

    extracted_neighbors = []
    for feat in record.features:
        if feat.type == "CDS":
            f_start = int(feat.location.start)
            f_end = int(feat.location.end)
            
            # Simple spatial exclusion to prevent self-overlap issues
            if not (a_start - window <= f_start <= a_end + window):
                continue
                
            n_strand = feat.location.strand
            direction = "Same" if n_strand == a_strand else "Oppose"
            
            # Signed Distance Calculation (Center of neighbor - Anchor position)
            dist = int(((f_start + f_end) / 2) - ((a_start + a_end) / 2))
            
            product = feat.qualifiers.get("product", ["unknown"])[0]
            translation = feat.qualifiers.get("translation", [""])[0]
            
            # Find Uniprot ID of neighbor if available
            n_pid = "unknown"
            for ref in feat.qualifiers.get("db_xref", []):
                if "UniProtKB/TrEMBL" in ref or "UniProtKB/Swiss-Prot" in ref:
                    n_pid = ref.split(":")[1]

            extracted_neighbors.append({
                "pid": n_pid,
                "product": product,
                "dist": dist,
                "direction": direction,
                "seq": translation
            })

    return pid, instance_data, extracted_neighbors

async def run_pipeline_async(homolog_ids):
    """
    Orchestrates the asynchronous fetching of UniProt and ENA data.
    """
    client = UniProtENAClient(max_concurrent=15)
    
    async with aiohttp.ClientSession() as session:
        tasks = [process_homolog(session, client, pid) for pid in homolog_ids]
        results = await tqdm.gather(*tasks, desc="Fetching Genomic Contexts")
        return results

def run_gnn_scout(input_fasta_path, hit_limit=500):
    discovery = DiscoveryEngine()
    db = GNNDB()

    with open(input_fasta_path, "r") as f:
        fasta_str = f.read()

    # 1. Initial Profile Search (now via HMMER/UniProt)
    homolog_ids, scores = discovery.find_homologs(fasta_str, hit_limit=hit_limit)

    if scores:
        display_search_metrics(scores)

    # 2. Asynchronous Context Fetching
    print(f"\nInitiating async retrieval for {len(homolog_ids)} targets...")
    results = asyncio.run(run_pipeline_async(homolog_ids))

    # 3. Synchronous Database Writes (Protects SQLite integrity)
    print("\nCommitting neighborhood data to database...")
    for pid, instance_data, neighbors in results:
        if not instance_data or not neighbors:
            continue
            
        instance_id = db.add_instance(
            pid, 
            instance_data["embl_acc"], 
            instance_data["start"], 
            instance_data["end"], 
            instance_data["strand"]
        )
        
        for n in neighbors:
            db.add_neighbor(
                instance_id, 
                n["pid"], 
                n["product"], 
                n["dist"], 
                n["direction"], 
                n["seq"]
            )

    print("Data extraction complete.")

    # 4. REPORT
    generate_summary_report()

def get_maximal_modules(raw_results, min_support=3, min_module_size=4):
    """
    Identifies functional modules and removes redundant subsets.
    """
    neighborhoods = defaultdict(set)
    for instance_id, product in raw_results:
        neighborhoods[instance_id].add(product)

    itemset_counts = Counter()
    
    for neighbors in neighborhoods.values():
        if len(neighbors) < min_module_size:
            continue
            
        sorted_neighbors = sorted(list(neighbors))
        if len(sorted_neighbors) > 20: 
             continue 

        for r in range(min_module_size, len(sorted_neighbors) + 1):
             for combo in combinations(sorted_neighbors, r):
                 itemset_counts[combo] += 1

    candidates = []
    for genes, count in itemset_counts.items():
        if count >= min_support:
            candidates.append({'genes': set(genes), 'count': count, 'display': genes})

    candidates.sort(key=lambda x: len(x['genes']), reverse=True)
    final_modules = []
    
    for candidate in candidates:
        is_redundant = False
        
        for larger_mod in final_modules:
            if candidate['genes'].issubset(larger_mod['genes']):
                if candidate['count'] == larger_mod['count']:
                    is_redundant = True
                    break
        
        if not is_redundant:
            final_modules.append(candidate)

    final_modules.sort(key=lambda x: (-x['count'], -len(x['genes'])))
    top_modules = final_modules[:10]
    
    return top_modules

def display_search_metrics(scores):
    """
    Displays distribution of search scores.
    Adapted from percentage identity to handle HMMER bit scores.
    """
    hi = max(scores)
    lo = min(scores)
    avg = sum(scores) / len(scores)

    print("\n" + "-"*60)
    print(f"Initial Search Info (n={len(scores)})")
    print(f"Score Range: {lo:.1f} - {hi:.1f} | Average Score: {avg:.1f}")
    print("-" * 60 + "\n")

def generate_summary_report():
    print("\n" + "="*145)
    print(f"{'GENOMIC NEIGHBORHOOD & CO-OCCURRENCE REPORT':^145}")
    print("="*145)
    
    conn = sqlite3.connect("data/GNN.db")
    cursor = conn.cursor()
    
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
    
    print(f"{'Neighbor Product (Example Accession)':<65} | {'Strand':<8} | {'Freq':<5} | {'Avg Dist':<12} | {'Spread':<8} | {'Len'}")
    print("-" * 145)
    
    for row in results:
        display_name = f"{row[0]} ({row[1]})"
        display_name = (display_name[:62] + '..') if len(display_name) > 65 else display_name
        print(f"{display_name:<65} | {row[2]:<8} | {row[3]:<5} | {row[4]:>9} bp | {row[5]:>6} bp | {row[6]}aa")

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

    raw_query = """
    SELECT instance_id, product
    FROM neighbors
    ORDER BY instance_id
    """
    cursor.execute(raw_query)
    raw_results = cursor.fetchall()

    modules = get_maximal_modules(raw_results, min_support=3, min_module_size=3)

    print("\n" + "="*80)
    print("TOP 10 VERIFIED FUNCTIONAL MODULES BY FREQ (Maximal Sets Only)")
    print("="*80)

    if not modules:
        print("No modules found matching criteria.")
    else:
        for i, mod in enumerate(modules, 1):
            print(f"{i}. [Freq: {mod['count']}] Size: {len(mod['genes'])}")
            clean_genes = ', '.join(sorted(list(mod['genes'])))
            print(f"   Genes: {clean_genes}")
            print("-" * 40)
            
    conn.close()
    print("="*145)

if __name__ == "__main__":
    # Ensure standard execution expects at least 1 sys.argv argument if run via CLI
    input_file = sys.argv[1] if len(sys.argv) > 1 else "input.fasta"
    run_gnn_scout(input_file, hit_limit=500)