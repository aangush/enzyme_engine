import os
from Bio import SeqIO
from collections import defaultdict

RAW_DIR = "data/raw"

def analyze_neighborhoods(anchor_accession_substring="hydrolase"):
    """
    Parses cached GenBank files to find common genomic neighbors.
    """
    # Dictionary to store neighbor stats: { "Product Name": [list of distances] }
    neighbor_data = defaultdict(list)
    total_genomes = 0

    # Loop through every cached GenBank file
    for filename in os.listdir(RAW_DIR):
        if not filename.endswith(".gbk"):
            continue
            
        path = os.path.join(RAW_DIR, filename)
        record = SeqIO.read(path, "genbank")
        total_genomes += 1
        
        # 1. Find the 'Anchor' gene in this specific slice
        anchor_feature = None
        cds_features = [f for f in record.features if f.type == "CDS"]
        
        for feat in cds_features:
            product = feat.qualifiers.get("product", [""])[0].lower()
            prot_id = feat.qualifiers.get("protein_id", [""])[0]
            
            # We identify the anchor by checking if it matches our query protein
            if anchor_accession_substring.lower() in product or anchor_accession_substring in prot_id:
                anchor_feature = feat
                break
        
        if not anchor_feature:
            continue

        anchor_center = (anchor_feature.location.start + anchor_feature.location.end) / 2

        # 2. Calculate distance to every other neighbor
        for feat in cds_features:
            if feat == anchor_feature:
                continue
                
            product = feat.qualifiers.get("product", ["hypothetical protein"])[0]
            feat_center = (feat.location.start + feat.location.end) / 2
            
            # Absolute distance in base pairs
            distance = abs(anchor_center - feat_center)
            
            neighbor_data[product].append(distance)

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
    data, count = analyze_neighborhoods("hydrolase")
    print_gnn_report(data, count)