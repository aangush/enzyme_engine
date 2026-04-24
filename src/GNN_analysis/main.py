import os
import sys
import sqlite3
import asyncio
import aiohttp
from collections import defaultdict, Counter
from itertools import combinations
try:
    from tqdm.asyncio import tqdm
except ImportError:
    print("Please install tqdm: pip install tqdm")
    sys.exit(1)

from discovery import DiscoveryEngine
from ena_client import UniProtENAClient
from db_manager import GNNDB
from domain_analyzer import DomainAnalyzer

async def process_homolog(session, client, pid, window_cds=15, required_pfams=None):
    """
    Async worker to fetch metadata, ENA context, and return structured neighbor data.
    """
    metadata, locations = await client.fetch_uniprot_data(session, pid)
    
    if not locations or not metadata:
        return pid, []

    results = []
    
    # Process up to 3 contigs per UniProt ID to resolve the statistical bottleneck
    for loc in locations[:3]:
        embl_acc = loc['embl_acc']
        record = await client.fetch_ena_neighborhood(session, embl_acc)
        if not record:
            continue

        # Isolate CDS features to normalize window size across genomes
        cds_features = [feat for feat in record.features if feat.type == "CDS"]

        anchor_idx = None
        for i, feat in enumerate(cds_features):
            uniprot_refs = [ref for ref in feat.qualifiers.get("db_xref", []) if "UniProtKB" in ref]
            if any(pid in ref for ref in uniprot_refs):
                anchor_idx = i
                break

        if anchor_idx is None:
            continue

        a_feat = cds_features[anchor_idx]
        a_start = int(a_feat.location.start)
        a_end = int(a_feat.location.end)
        a_strand = a_feat.location.strand

        instance_data = {
            "embl_acc": embl_acc,
            "start": a_start,
            "end": a_end,
            "strand": a_strand
        }

        extracted_neighbors = []
        locus_pfams = set()
        
        start_idx = max(0, anchor_idx - window_cds)
        end_idx = min(len(cds_features), anchor_idx + window_cds + 1)

        for i in range(start_idx, end_idx):
            if i == anchor_idx:
                continue
                
            feat = cds_features[i]
            f_start = int(feat.location.start)
            f_end = int(feat.location.end)
            
            n_strand = feat.location.strand
            direction = "Same" if n_strand == a_strand else "Oppose"
            dist = int(((f_start + f_end) / 2) - ((a_start + a_end) / 2))
            
            product = feat.qualifiers.get("product", ["unknown"])[0]
            translation = feat.qualifiers.get("translation", [""])[0]
            
            # Standardize nodes using ENA database cross-references
            pfam_ids = [ref.split(":")[1] for ref in feat.qualifiers.get("db_xref", []) if "Pfam" in ref or "PANTHER" in ref or "InterPro" in ref]
            locus_pfams.update(pfam_ids)
            pfam_str = ",".join(pfam_ids)
            
            n_pid = "unknown"
            for ref in feat.qualifiers.get("db_xref", []):
                if "UniProtKB/TrEMBL" in ref or "UniProtKB/Swiss-Prot" in ref:
                    n_pid = ref.split(":")[1]

            # Proxy filter for eukaryotic/membrane compatibility
            notes = " ".join(feat.qualifiers.get("note", [])).lower()
            is_membrane = 1 if ("membrane" in notes or "transporter" in notes or "membrane" in product.lower()) else 0

            extracted_neighbors.append({
                "pid": n_pid,
                "product": product,
                "dist": dist,
                "direction": direction,
                "seq": translation,
                "pfam_ids": pfam_str,
                "is_membrane": is_membrane
            })

        # Multi-anchor locus requirement filter
        if required_pfams:
            overlap = locus_pfams.intersection(set(required_pfams))
            if len(overlap) < len(required_pfams):
                continue # Discard this contig; it lacks the necessary pathway components

        results.append((instance_data, extracted_neighbors))

    return pid, results

async def run_pipeline_async(homolog_ids, required_pfams=None):
    client = UniProtENAClient(max_concurrent=15)
    async with aiohttp.ClientSession() as session:
        tasks = [process_homolog(session, client, pid, required_pfams=required_pfams) for pid in homolog_ids]
        results = await tqdm.gather(*tasks, desc="Fetching Genomic Contexts")
        return results

def run_gnn_scout(input_fasta_path, hit_limit=500, required_pfams=None):
    discovery = DiscoveryEngine()
    db = GNNDB()

    with open(input_fasta_path, "r") as f:
        fasta_str = f.read()

    homolog_ids, scores = discovery.find_homologs(fasta_str, hit_limit=hit_limit)

    print(f"\nInitiating async retrieval for {len(homolog_ids)} targets...")
    # Pass required pfams to async runner
    results = asyncio.run(run_pipeline_async(homolog_ids, required_pfams=required_pfams))

    print("\nCommitting neighborhood data to database...")
    for pid, locus_results in results:
        for instance_data, neighbors in locus_results:
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
                    n["seq"],
                    n["pfam_ids"],
                    n["is_membrane"]
                )

    print("Data extraction complete.")
    analyzer = DomainAnalyzer(db_path="data/GNN.db")
    analyzer.analyze_hypotheticals()
    generate_summary_report()

def get_maximal_modules(raw_results, min_support=3, min_module_size=4):
    neighborhoods = defaultdict(set)
    for instance_id, identifier in raw_results:
        # Operating on the standardized identifier (Pfam or product string)
        neighborhoods[instance_id].add(identifier)

    itemset_counts = Counter()
    for neighbors in neighborhoods.values():
        if len(neighbors) < min_module_size: continue
        sorted_neighbors = sorted(list(neighbors))
        if len(sorted_neighbors) > 25: continue 

        for r in range(min_module_size, len(sorted_neighbors) + 1):
             for combo in combinations(sorted_neighbors, r):
                 itemset_counts[combo] += 1

    candidates = [{'genes': set(genes), 'count': count, 'display': genes} 
                  for genes, count in itemset_counts.items() if count >= min_support]

    candidates.sort(key=lambda x: len(x['genes']), reverse=True)
    final_modules = []
    
    for candidate in candidates:
        is_redundant = False
        for larger_mod in final_modules:
            if candidate['genes'].issubset(larger_mod['genes']) and candidate['count'] == larger_mod['count']:
                is_redundant = True
                break
        if not is_redundant:
            final_modules.append(candidate)

    final_modules.sort(key=lambda x: (-x['count'], -len(x['genes'])))
    return final_modules[:10]

def generate_summary_report():
    print("\n" + "="*80)
    print("TOP 10 VERIFIED FUNCTIONAL MODULES BY FREQ (Maximal Sets Only)")
    print("="*80)
    
    conn = sqlite3.connect("data/GNN.db")
    cursor = conn.cursor()

    # Standardization mapping: Use Pfam ID as node identity if available, else product name
    raw_query = """
    SELECT instance_id, 
           CASE WHEN pfam_ids != '' THEN pfam_ids ELSE product END as identifier
    FROM neighbors
    ORDER BY instance_id
    """
    cursor.execute(raw_query)
    raw_results = cursor.fetchall()

    modules = get_maximal_modules(raw_results, min_support=3, min_module_size=3)

    if not modules:
        print("No modules found matching criteria.")
    else:
        for i, mod in enumerate(modules, 1):
            print(f"{i}. [Freq: {mod['count']}] Size: {len(mod['genes'])}")
            clean_genes = ', '.join(sorted(list(mod['genes'])))
            print(f"   Network Nodes: {clean_genes}")
            print("-" * 40)
            
    conn.close()

if __name__ == "__main__":
    input_file = sys.argv[1] if len(sys.argv) > 1 else "input.fasta"
    
    # Example usage for the multi-anchor constraint:
    # Set this list to the known Pfam domains of your 8-enzyme pathway.
    # The pipeline will only return genomic windows containing these domains.
    mp_pathway_pfams = ["IPR001559", "IPR050631", "IPR037513", "IPR050770", "PTHR11699"]
    
    run_gnn_scout(input_file, hit_limit=500, required_pfams=mp_pathway_pfams)