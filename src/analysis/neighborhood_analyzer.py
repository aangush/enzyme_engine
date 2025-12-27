import os
from Bio import SeqIO
from collections import defaultdict

RAW_DIR = "data/raw"

def analyze_neighborhoods(center_point=10000, min_length=15000):
    """
    Parses cached GenBank files. 
    Identifies anchor by position (center) and treats hypothetical proteins as individuals.
    """
    neighbor_data = defaultdict(list)
    total_genomes = 0

    if not os.path.exists(RAW_DIR):
        print(f"Error: Directory {RAW_DIR} not found.")
        return neighbor_data, 0

    for filename in os.listdir(RAW_DIR):
        if not filename.endswith(".gbk"):
            continue
            
        path = os.path.join(RAW_DIR, filename)
        
        # FIX: Load the record FIRST, then you can check its length
        try:
            record = SeqIO.read(path, "genbank")
        except Exception as e:
            print(f"Skipping {filename}: Could not read file. {e}")
            continue
        
        # Now we can safely check the length
        if len(record.seq) < min_length:
            print(f"Skipping {filename}: Sequence too short ({len(record.seq)} bp)")
            continue
        
        total_genomes += 1
        cds_features = [f for f in record.features if f.type == "CDS"]
        
        # 1. FIND ANCHOR BY POSITION
        anchor_feature = None
        for feat in cds_features:
            if feat.location.start <= center_point <= feat.location.end:
                anchor_feature = feat
                break
        
        # Fallback if center point isn't exactly inside a gene
        if not anchor_feature and cds_features:
            cds_features.sort(key=lambda f: abs(((f.location.start + f.location.end)/2) - center_point))
            anchor_feature = cds_features[0]

        if not anchor_feature:
            continue

        # --- GET ANCHOR STRAND AND CENTER ---
        anchor_strand = anchor_feature.location.strand
        anchor_center = (anchor_feature.location.start + anchor_feature.location.end) / 2

        # 2. PROCESS NEIGHBORS
        for feat in cds_features:
            if feat == anchor_feature:
                continue

            # Clean up the product name
            raw_product = feat.qualifiers.get("product", ["hypothetical protein"])[0]
            product = raw_product.strip().lower()
            
            # Identify label
            if "hypothetical" in product:
                prot_id = feat.qualifiers.get("protein_id", ["no_id"])[0]
                label = f"Hypothetical ({prot_id})"
            else:
                label = product.capitalize()

            # --- STRAND LOGIC ---
            neighbor_strand = feat.location.strand
            direction = "Same" if anchor_strand == neighbor_strand else "Oppose"
            final_label = f"{label} [{direction}]"    

            feat_center = (feat.location.start + feat.location.end) / 2
            distance = abs(anchor_center - feat_center)
            
            neighbor_data[final_label].append(distance)

    return neighbor_data, total_genomes

def print_gnn_report(neighbor_data, total_genomes):
    if total_genomes == 0:
        print("No valid neighborhoods found to analyze.")
        return

    print(f"\n{'='*80}")
    print(f"GENOMIC NEIGHBORHOOD NETWORK REPORT (n={total_genomes} genomes)")
    print(f"{'='*80}")
    print(f"{'Neighboring Product':<65} {'Freq':<10} {'Avg Dist (bp)':<15}")
    print(f"{'-'*80}")

    # Sort by frequency
    sorted_neighbors = sorted(neighbor_data.items(), 
                              key=lambda x: len(x[1]), 
                              reverse=True)

    for product, distances in sorted_neighbors[:20]: 
        freq = len(distances)
        avg_dist = sum(distances) / freq
        # We only care about neighbors that appear in more than 1 genome
        if freq > 1:
            print(f"{product[:62]:<65} {freq:<10} {int(avg_dist):<15}")

if __name__ == "__main__":
    data, count = analyze_neighborhoods()
    print_gnn_report(data, count)