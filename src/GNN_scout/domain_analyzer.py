import sqlite3
import subprocess
import os

class DomainAnalyzer:
    def __init__(self, db_path="data/scout.db", pfam_path="data/pfam/Pfam-A.hmm"):
        self.db_path = db_path
        self.pfam_path = pfam_path

    def analyze_hypotheticals(self):
        print("Starting Hidden Domain Analysis...")
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        # 1. Get unique hypothetical sequences
        cursor.execute("SELECT DISTINCT protein_id, sequence FROM neighbors WHERE product LIKE '%hypothetical%' AND sequence != ''")
        hypotheticals = cursor.fetchall()
        
        print(f"Found {len(hypotheticals)} unique hypothetical proteins to analyze.")

        for prot_id, seq in hypotheticals:
            domain = self._scan_sequence(prot_id, seq)
            if domain:
                print(f"Identified: {prot_id} -> {domain}")
                # Update the product name in the DB to include the domain
                new_product = f"Hypothetical ({domain})"
                cursor.execute("UPDATE neighbors SET product = ? WHERE protein_id = ?", (new_product, prot_id))
        
        conn.commit()
        conn.close()
        print("Domain analysis complete.")

    def _scan_sequence(self, prot_id, sequence):
        # Create a temporary FASTA file
        fasta_path = f"data/temp_{prot_id}.fasta"
        with open(fasta_path, "w") as f:
            f.write(f">{prot_id}\n{sequence}\n")

        # Run hmmscan (local)
        try:
            cmd = ["hmmscan", "--tblout", f"data/temp_{prot_id}.txt", "--noali", "-E", "1e-5", self.pfam_path, fasta_path]
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            # Parse the tblout file
            with open(f"data/temp_{prot_id}.txt", "r") as res:
                for line in res:
                    if not line.startswith("#"):
                        parts = line.split()
                        domain_name = parts[0] # The Pfam entry name
                        return domain_name
        except Exception as e:
            print(f"Error scanning {prot_id}: {e}")
        finally:
            # Clean up temp files
            if os.path.exists(fasta_path): os.remove(fasta_path)
            if os.path.exists(f"data/temp_{prot_id}.txt"): os.remove(f"data/temp_{prot_id}.txt")
        
        return None