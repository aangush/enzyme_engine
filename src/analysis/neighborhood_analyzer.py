import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

# Setup directories
RAW_DIR = "data/raw"
INTERIM_DIR = "data/interim"
os.makedirs(INTERIM_DIR, exist_ok=True)

def sequence_identity(seq1, seq2):
    """Calculates basic sequence identity to group similar hypothetical proteins."""
    if not seq1 or not seq2: return 0
    shorter = min(len(seq1), len(seq2))
    matches = sum(1 for i in range(shorter) if seq1[i] == seq2[i])
    return (matches / max(len(seq1), len(seq2))) * 100

def cluster_hypotheticals(hypo_data, threshold=90):
    """Groups protein IDs that share >90% sequence identity."""
    clusters = {}
    cluster_count = 0
    ids = list(hypo_data.keys())
    
    for i in range(len(ids)):
        if ids[i] in clusters: continue
        
        cluster_count += 1
        current_cluster_id = f"HypoCluster_{cluster_count:02}"
        clusters[ids[i]] = current_cluster_id
        
        for j in range(i + 1, len(ids)):
            if ids[j] in clusters: continue
            identity = sequence_identity(hypo_data[ids[i]], hypo_data[ids[j]])
            if identity >= threshold:
                clusters[ids[j]] = current_cluster_id
                
    return clusters

def analyze_neighborhoods(center_point=10000, min_length=15000):
    """
    Parses cached GenBank files, identifies anchors by position, 
    and clusters hypothetical proteins by sequence similarity.
    """
    neighbor_data = defaultdict(list)
    hypo_sequences = {} # {prot_id: sequence}
    raw_observations = [] # Temp storage for hypothetical hits
    total_genomes = 0

    if not os.path.exists(RAW_DIR):
        print(f"Error: Directory {RAW_DIR} not found.")
        return neighbor_data, 0, []

    for filename in os.listdir(RAW_DIR):
        if not filename.endswith(".gbk"): continue
        path = os.path.join(RAW_DIR, filename)
        
        try:
            record = SeqIO.read(path, "genbank")
        except Exception: continue
        
        if len(record.seq) < min_length:
            continue
        
        total_genomes += 1
        cds_features = [f for f in record.features if f.type == "CDS"]
        
        # 1. Find Anchor by Position
        anchor_feature = next((f for f in cds_features if f.location.start <= center_point <= f.location.end), None)
        if not anchor_feature and cds_features:
            cds_features.sort(key=lambda f: abs(((f.location.start + f.location.end)/2) - center_point))
            anchor_feature = cds_features[0]

        if not anchor_feature: continue
        anchor_strand = anchor_feature.location.strand
        anchor_center = (anchor_feature.location.start + anchor_feature.location.end) / 2

        # 2. Process Neighbors
        for feat in cds_features:
            if feat == anchor_feature: continue

            raw_product = feat.qualifiers.get("product", ["hypothetical protein"])[0]
            product = raw_product.strip().lower()
            prot_id = feat.qualifiers.get("protein_id", ["no_id"])[0]
            
            direction = "Same" if anchor_strand == feat.location.strand else "Oppose"
            dist = abs(anchor_center - ((feat.location.start + feat.location.end) / 2))

            if "hypothetical" in product:
                translation = feat.qualifiers.get("translation", [""])[0]
                if translation:
                    hypo_sequences[prot_id] = translation
                # Store for later clustering
                raw_observations.append({"id": prot_id, "dist": dist, "dir": direction})
            else:
                label = f"{product.capitalize()} [{direction}]"
                neighbor_data[label].append(dist)

    # --- CLUSTERING STEP ---
    hypo_clusters = cluster_hypotheticals(hypo_sequences)
    
    # Add clustered hypotheticals back into the neighbor_data
    for obs in raw_observations:
        c_id = hypo_clusters.get(obs["id"], "Unique_Hypo")
        label = f"{c_id} [{obs['dir']}]"
        neighbor_data[label].append(obs["dist"])

    # Prepare export records
    seq_records = [SeqRecord(Seq(s), id=pid, description=f"Cluster:{hypo_clusters.get(pid)}") 
                   for pid, s in hypo_sequences.items()]

    return neighbor_data, total_genomes, seq_records

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

def export_hypotheticals(seq_records):
    if not seq_records: return
    fasta_path = os.path.join(INTERIM_DIR, "neighborhood_hypotheticals.fasta")
    SeqIO.write(seq_records, fasta_path, "fasta")
    print(f"🧬 Exported {len(seq_records)} hypothetical sequences to {fasta_path}")

if __name__ == "__main__":
    data, count, hypos = analyze_neighborhoods()
    print_gnn_report(data, count)
    export_hypotheticals(hypos)