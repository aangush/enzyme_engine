import sqlite3
import subprocess
import os
from difflib import SequenceMatcher

class DomainAnalyzer:
    def __init__(self, db_path="data/scout.db", pfam_path="data/pfam/Pfam-A.hmm"):
        self.db_path = db_path
        self.pfam_path = pfam_path

    def _get_identity(self, a, b):
        """Calculates sequence identity ratio."""
        return SequenceMatcher(None, a, b).ratio()

    def analyze_hypotheticals(self):
        print("\n" + "="*70)
        print("SEQUENCE-BASED NO PFAM CLUSTERING (70% IDENTITY)")
        print("="*70)
        
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        # Fetch targets: Must have sequence, must not have been processed
        cursor.execute("""
            SELECT id, protein_id, sequence 
            FROM neighbors 
            WHERE (LOWER(product) LIKE '%hypothetical%' OR LOWER(product) = 'unknown') 
            AND sequence != '' AND product NOT LIKE '%(%'
        """)
        targets = cursor.fetchall()
        
        if not targets:
            print("No new un-analyzed sequences found.")
            conn.close()
            return

        # key: cluster_index -> representative_sequence
        mystery_representatives = {} 
        next_cluster_id = 1

        print(f"Processing {len(targets)} targets...")

        for db_id, prot_id, seq in targets:
            domain = self._scan_sequence(prot_id, seq)
            
            if domain:
                new_product = f"Hypothetical ({domain})"
            else:
                # SEQUENCE IDENTITY CLUSTERING
                assigned_cluster = None
                
                for cid, rep_seq in mystery_representatives.items():
                    if self._get_identity(seq, rep_seq) >= 0.70:
                        assigned_cluster = cid
                        break
                
                if assigned_cluster is None:
                    assigned_cluster = next_cluster_id
                    mystery_representatives[assigned_cluster] = seq
                    next_cluster_id += 1
                
                # Tagging by Cluster ID instead of distance
                new_product = f"Hypothetical (Mystery Cluster {assigned_cluster})"
            
            cursor.execute("UPDATE neighbors SET product = ? WHERE id = ?", (new_product, db_id))
        
        conn.commit()
        conn.close()
        print(f"Complete. Grouped mystery proteins into {next_cluster_id - 1} sequence clusters.")

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