import sqlite3
import subprocess
import os
from difflib import SequenceMatcher

class DomainAnalyzer:
    def __init__(self, db_path="data/scout.db", pfam_path="data/pfam/Pfam-A.hmm"):
        self.db_path = db_path
        self.pfam_path = pfam_path

    def _get_identity(self, a, b):
        """Quick sequence identity ratio."""
        return SequenceMatcher(None, a, b).ratio()

    def analyze_hypotheticals(self):
        print("\n" + "="*70)
        print("HYBRID SYNTENY-IDENTITY ANALYSIS")
        print("="*70)
        
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        # 1. Fetch un-analyzed targets
        cursor.execute("""
            SELECT id, protein_id, sequence, distance_bp 
            FROM neighbors 
            WHERE (LOWER(product) LIKE '%hypothetical%' OR LOWER(product) = 'unknown') 
            AND sequence != '' AND product NOT LIKE '%(%'
        """)
        targets = cursor.fetchall()
        
        if not targets:
            print("No new sequences found.")
            conn.close()
            return

        # Tracking for identity-based clustering of mysteries
        # Key: (bin_dist, len_bucket) -> List of sequences
        mystery_clusters = {} 

        for db_id, prot_id, seq, dist in targets:
            domain = self._scan_sequence(prot_id, seq)
            
            if domain:
                new_product = f"Hypothetical ({domain})"
            else:
                # SPATIAL + LENGTH BINNING
                binned_dist = round(dist / 1000) * 1000
                len_bucket = round(len(seq) / 50) * 50
                cluster_key = (binned_dist, len_bucket)
                
                # IDENTITY CHECK
                cluster_match = False
                if cluster_key in mystery_clusters:
                    for rep_seq in mystery_clusters[cluster_key]:
                        if self._get_identity(seq, rep_seq) > 0.90:
                            cluster_match = True
                            break
                
                if not cluster_match:
                    if cluster_key not in mystery_clusters:
                        mystery_clusters[cluster_key] = []
                    mystery_clusters[cluster_key].append(seq)
                
                prefix = "+" if binned_dist >= 0 else ""
                new_product = f"Hypothetical (No Pfam Hit) @ {prefix}{binned_dist}bp | ~{len_bucket}aa"
            
            cursor.execute("UPDATE neighbors SET product = ? WHERE id = ?", (new_product, db_id))
        
        conn.commit()
        conn.close()
        print("Analysis complete. Mystery clusters merged by location and identity.")

    def _scan_sequence(self, prot_id, sequence):
        fasta_path = f"data/temp_{prot_id}.fasta"
        out_path = f"data/temp_{prot_id}.txt"
        with open(fasta_path, "w") as f: f.write(f">{prot_id}\n{sequence}\n")
        try:
            cmd = ["hmmscan", "--tblout", out_path, "--noali", "-E", "1e-5", self.pfam_path, fasta_path]
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if os.path.exists(out_path):
                with open(out_path, "r") as res:
                    for line in res:
                        if not line.startswith("#"):
                            return line.split()[0]
        except Exception: return None
        finally:
            for p in [fasta_path, out_path]:
                if os.path.exists(p): os.remove(p)
        return None