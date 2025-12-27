import os
from Bio import SeqIO
from collections import defaultdict

RAW_DIR = "data/raw"

def analyze_neighborhoods(center_point=10000):
    """
    Parses cached GenBank files. 
    Identifies anchor by position (center) and treats hypothetical proteins as individuals.
    """
    neighbor_data = defaultdict(list)
    total_genomes = 0

    for filename in os.listdir(RAW_DIR):
        if not filename.endswith(".gbk"):
            continue
            
        path = os.path.join(RAW_DIR, filename)
        record = SeqIO.read(path, "genbank")
        total_genomes += 1
        
        cds_features = [f for f in record.features if f.type == "CDS"]
        
        # 1. FIND ANCHOR BY POSITION
        # The anchor is the gene that contains the center point (usually 10,000)
        anchor_feature = None
        for feat in cds_features:
            if feat.location.start <= center_point <= feat.location.end:
                anchor_feature = feat
                break
        
        if not anchor_feature:
            # Fallback: find the one closest to the center
            cds_features.sort(key=lambda f: abs(((f.location.start + f.location.end)/2) - center_point))
            anchor_feature = cds_features[0]

        anchor_center = (anchor_feature.location.start + anchor_feature.location.end) / 2

        # 2. PROCESS NEIGHBORS
        for feat in cds_features:
            if feat == anchor_feature:
                continue
                
            # Clean up the product name
            raw_product = feat.qualifiers.get("product", ["hypothetical protein"])[0]
            product = raw_product.strip().lower() # Standardize case
            
            # If it's hypothetical, attach the ID so they aren't all 'blobbed' together
            if "hypothetical" in product:
                prot_id = feat.qualifiers.get("protein_id", ["no_id"])[0]
                label = f"Hypothetical ({prot_id})"
            else:
                label = product.capitalize()

            feat_center = (feat.location.start + feat.location.end) / 2
            distance = abs(anchor_center - feat_center)
            
            neighbor_data[label].append(distance)

    return neighbor_data, total_genomes

def print_gnn_report(neighbor_data, total_genomes):
    print(f"\n{'='*80}")
    print(f"GENOMIC NEIGHBORHOOD NETWORK REPORT (n={total_genomes} genomes)")
    print(f"{'='*80}")
    print(f"{'Neighboring Product':<45} {'Freq':<10} {'Avg Dist (bp)':<15}")
    print(f"{'-'*80}")

    # Sort by frequency (most common neighbors first)
    sorted_neighbors = sorted(neighbor_data.items(), 
                              key=lambda x: len(x[1]), 
                              reverse=True)

    for product, distances in sorted_neighbors[:15]: # Top 15 neighbors
        freq = len(distances)
        avg_dist = sum(distances) / freq
        # We only care about neighbors that appear in more than 1 genome
        if freq > 1:
            print(f"{product[:43]:<45} {freq:<10} {int(avg_dist):<15}")

if __name__ == "__main__":
    # Run the analysis
    data, count = analyze_neighborhoods()
    print_gnn_report(data, count)