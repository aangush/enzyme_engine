import sqlite3
import subprocess
import os

class DomainAnalyzer:
    def __init__(self, db_path="data/scout.db", pfam_path="data/pfam/Pfam-A.hmm"):
        self.db_path = db_path
        self.pfam_path = pfam_path

    def analyze_hypotheticals(self):
        """
        Identifies unannotated proteins and attempts to assign functional domains
        using local HMMER scans against the Pfam database.
        """
        print("\n" + "="*60)
        print("STARTING HIDDEN DOMAIN ANALYSIS")
        print("="*60)
        
        if not os.path.exists(self.pfam_path):
            print(f"Error: Pfam database not found at {self.pfam_path}")
            print("Please ensure you have run 'hmmpress' on the Pfam-A.hmm file.")
            return

        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        # Target: Products that are some variation of 'hypothetical' or 'unknown'
        # Criteria: Must have a sequence, and must NOT have been analyzed yet (no parentheses)
        cursor.execute("""
            SELECT DISTINCT protein_id, sequence 
            FROM neighbors 
            WHERE (LOWER(product) LIKE '%hypothetical%' OR LOWER(product) = 'unknown') 
            AND sequence != '' 
            AND product NOT LIKE '%(%'
        """)
        
        targets = cursor.fetchall()
        total = len(targets)
        
        if total == 0:
            print("No un-analyzed hypothetical proteins found in database.")
            conn.close()
            return

        print(f"Found {total} unique sequences to analyze. Processing...")

        hits = 0
        for idx, (prot_id, seq) in enumerate(targets, 1):
            # Simple progress tracker
            if idx % 5 == 0 or idx == total:
                print(f"Progress: {idx}/{total} sequences scanned...")

            domain = self._scan_sequence(prot_id, seq)
            
            if domain:
                new_product = f"Hypothetical ({domain})"
                hits += 1
            else:
                # Labeling it ensures we don't waste time re-scanning it in future runs
                new_product = "Hypothetical (No Pfam Hit)"
            
            cursor.execute("UPDATE neighbors SET product = ? WHERE protein_id = ?", (new_product, prot_id))
        
        conn.commit()
        conn.close()
        print("="*60)
        print(f"ANALYSIS COMPLETE: Identified {hits} new domains across {total} targets.")
        print("="*60)

    def _scan_sequence(self, prot_id, sequence):
        """Runs hmmscan for a single sequence and returns the top hit domain name."""
        fasta_path = f"data/temp_{prot_id}.fasta"
        out_path = f"data/temp_{prot_id}.txt"
        
        # 1. Create temporary FASTA
        with open(fasta_path, "w") as f:
            f.write(f">{prot_id}\n{sequence}\n")

        try:
            # 2. Run hmmscan
            # -E 1e-5: Sets a scientifically robust E-value cutoff
            # --domtblout: Is excellent for domain-based parsing
            cmd = [
                "hmmscan", 
                "--tblout", out_path, 
                "--noali", 
                "-E", "1e-5", 
                self.pfam_path, 
                fasta_path
            ]
            
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            # 3. Parse result
            if os.path.exists(out_path):
                with open(out_path, "r") as res:
                    for line in res:
                        if not line.startswith("#"):
                            parts = line.split()
                            if len(parts) > 0:
                                return parts[0] # The name of the Pfam domain (e.g., 'Pirin')
                                
        except Exception as e:
            print(f"Error scanning {prot_id}: {e}")
        finally:
            # 4. Cleanup temporary files immediately
            self._cleanup(prot_id)
        
        return None

    def _cleanup(self, prot_id):
        """Removes temporary files to keep the data directory clean."""
        for ext in [".fasta", ".txt"]:
            path = f"data/temp_{prot_id}{ext}"
            if os.path.exists(path):
                os.remove(path)