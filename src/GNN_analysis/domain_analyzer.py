import sqlite3
import subprocess
import os
from Bio import Align

class DomainAnalyzer:
    def __init__(self, db_path="data/GNN.db", pfam_path="data/pfam/Pfam-A.hmm"):
        self.db_path = db_path
        self.pfam_path = pfam_path

        # Initialize aligner and matrix
        self.aligner = Align.PairwiseAligner(mode="global")
        self.aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")

    def _get_identity(self, a, b):
        raw_score = self.aligner.score(a, b)
        shorter = a if len(a) < len(b) else b
        max_score = self.aligner.score(shorter, shorter)
        if max_score == 0: return 0.0
        return raw_score / max_score

    def analyze_hypotheticals(self):
        print("\n" + "="*70)
        print("SEQUENCE & DOMAIN CLUSTERING FOR UNCHARACTERIZED PROTEINS")
        print("="*70)
        
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        # Added 'uncharacterized' to catch ENA standard nomenclature
        cursor.execute("""
            SELECT id, protein_id, sequence 
            FROM neighbors 
            WHERE (LOWER(product) LIKE '%hypothetical%' 
               OR LOWER(product) LIKE '%unknown%' 
               OR LOWER(product) LIKE '%uncharacterized%') 
            AND sequence != '' 
            AND product NOT LIKE '%(%'
        """)
        targets = cursor.fetchall()
        
        if not targets:
            print("No un-analyzed sequences found.")
            conn.close()
            return

        mystery_representatives = {} 
        next_cluster_id = 1

        print(f"Processing {len(targets)} uncharacterized targets via local HMMER/Pfam...")

        for db_id, prot_id, seq in targets:
            # Pass db_id to ensure unique temporary files
            domain = self._scan_sequence(db_id, seq)
            
            if domain:
                new_product = f"Hypothetical ({domain})"
            else:
                assigned_cluster = None
                
                for cid, rep_seq in mystery_representatives.items():
                    if self._get_identity(seq, rep_seq) >= 0.70:
                        assigned_cluster = cid
                        break
                
                if assigned_cluster is None:
                    assigned_cluster = next_cluster_id
                    mystery_representatives[assigned_cluster] = seq
                    next_cluster_id += 1
                
                new_product = f"Hypothetical (Mystery Cluster {assigned_cluster})"
            
            cursor.execute("UPDATE neighbors SET product = ? WHERE id = ?", (new_product, db_id))
        
        conn.commit()
        conn.close()
        print(f"Complete. Grouped mystery proteins into {next_cluster_id - 1} novel sequence clusters.")

    def _scan_sequence(self, db_id, sequence):
        # Use db_id to prevent file collisions when protein_id is 'unknown'
        fasta_path = f"data/temp_{db_id}.fasta"
        out_path = f"data/temp_{db_id}.txt"
        
        with open(fasta_path, "w") as f: 
            f.write(f">seq_{db_id}\n{sequence}\n")
            
        try:
            cmd = ["hmmscan", "--tblout", out_path, "--noali", "-E", "1e-5", self.pfam_path, fasta_path]
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if os.path.exists(out_path):
                with open(out_path, "r") as res:
                    for line in res:
                        if not line.startswith("#"):
                            return line.split()[0]
        except Exception: 
            return None
        finally:
            for p in [fasta_path, out_path]:
                if os.path.exists(p): os.remove(p)
        return None