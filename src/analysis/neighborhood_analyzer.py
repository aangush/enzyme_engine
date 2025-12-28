import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

RAW_DIR = "data/raw"
INTERIM_DIR = "data/interim"

# Ensure directories exist
os.makedirs(INTERIM_DIR, exist_ok=True)

def analyze_neighborhoods(center_point=10000, min_length=15000):
    neighbor_data = defaultdict(list)
    hypothetical_seqs = [] # List of SeqRecord objects for FASTA export
    total_genomes = 0

    if not os.path.exists(RAW_DIR):
        print(f"Error: Directory {RAW_DIR} not found.")
        return neighbor_data, 0, []

    for filename in os.listdir(RAW_DIR):
        if not filename.endswith(".gbk"):
            continue
            
        path = os.path.join(RAW_DIR, filename)
        
        try:
            record = SeqIO.read(path, "genbank")
        except Exception as e:
            print(f"Skipping {filename}: Could not read file. {e}")
            continue
        
        if len(record.seq) < min_length:
            continue
        
        total_genomes += 1
        cds_features = [f for f in record.features if f.type == "CDS"]
        
        # 1. FIND ANCHOR BY POSITION
        anchor_feature = None
        for feat in cds_features:
            if feat.location.start <= center_point <= feat.location.end:
                anchor_feature = feat
                break
        
        if not anchor_feature and cds_features:
            cds_features.sort(key=lambda f: abs(((f.location.start + f.location.end)/2) - center_point))
            anchor_feature = cds_features[0]

        if not anchor_feature:
            continue

        anchor_strand = anchor_feature.location.strand
        anchor_center = (anchor_feature.location.start + anchor_feature.location.end) / 2

        # 2. PROCESS NEIGHBORS
        for feat in cds_features:
            if feat == anchor_feature:
                continue

            raw_product = feat.qualifiers.get("product", ["hypothetical protein"])[0]
            product = raw_product.strip().lower()
            prot_id = feat.qualifiers.get("protein_id", ["no_id"])[0]
            
            # Labeling: We group all "hypothetical protein" together in the report
            # but we still save their unique sequences for your Diversity Engine
            if "hypothetical" in product:
                label = "Hypothetical protein"
                translation = feat.qualifiers.get("translation", [""])[0]
                if translation:
                    # Create a record for the FASTA file
                    hypothetical_seqs.append(SeqRecord(
                        Seq(translation), 
                        id=prot_id, 
                        description=f"from_{record.id}_neighbor_of_anchor"
                    ))
            else:
                label = product.capitalize()

            # Strand/Direction Logic
            neighbor_strand = feat.location.strand
            direction = "Same" if anchor_strand == neighbor_strand else "Oppose"
            final_label = f"{label} [{direction}]"    

            feat_center = (feat.location.start + feat.location.end) / 2
            distance = abs(anchor_center - feat_center)
            
            neighbor_data[final_label].append(distance)

    return neighbor_data, total_genomes, hypothetical_seqs

def export_hypotheticals(seq_records):
    """Saves hypothetical sequences to a FASTA file."""
    if not seq_records:
        return
    
    fasta_path = os.path.join(INTERIM_DIR, "neighborhood_hypotheticals.fasta")
    SeqIO.write(seq_records, fasta_path, "fasta")
    print(f"🧬 Exported {len(seq_records)} hypothetical sequences to {fasta_path}")

def print_gnn_report(neighbor_data, total_genomes):
    if total_genomes == 0:
        print("No valid neighborhoods found.")
        return

    name_w, freq_w, dist_w = 65, 10, 15
    total_w = name_w + freq_w + dist_w

    print(f"\n{'='*total_w}")
    print(f"GENOMIC NEIGHBORHOOD NETWORK REPORT (n={total_genomes} genomes)")
    print(f"{'='*total_w}")
    print(f"{'Neighboring Product':<{name_w}} {'Freq':<{freq_w}} {'Avg Dist (bp)':<{dist_w}}")
    print(f"{'-'*total_w}")

    sorted_neighbors = sorted(neighbor_data.items(), key=lambda x: len(x[1]), reverse=True)

    for product, distances in sorted_neighbors[:25]: 
        freq = len(distances)
        avg_dist = sum(distances) / freq
        if freq > 1:
            print(f"{product[:62]:<{name_w}} {freq:<{freq_w}} {int(avg_dist):<{dist_w}}")

if __name__ == "__main__":
    data, count, hypos = analyze_neighborhoods()
    print_gnn_report(data, count)
    export_hypotheticals(hypos)