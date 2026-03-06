import requests
import time
import sys

class DiscoveryEngine:
    def find_homologs(self, sequence_fasta, hit_limit=500):
        print(f"Submitting sequence search to EBI Job Dispatcher (NCBI BLAST+ against UniProtKB, limit={hit_limit})...")
        
        seq_clean = sequence_fasta.strip()
        if not seq_clean.startswith(">"):
            seq_clean = f">Query_Sequence\n{seq_clean}"

        submit_url = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run"
        payload = {
            "email": "enzyme_engine_user@gmail.com", 
            "program": "blastp",
            "stype": "protein",
            "sequence": seq_clean,
            "database": "uniprotkb",
            "alignments": str(hit_limit),
            "scores": str(hit_limit)
        }
        
        try:
            response = requests.post(submit_url, data=payload)
            if response.status_code != 200:
                print(f"API Error: Server rejected payload. Details: {response.text}")
                return [], []
                
            job_id = response.text.strip()
            print(f"Job submitted successfully. Job ID: {job_id}")
            print("Polling for results (this usually takes 1-3 minutes)...")
            
            status_url = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/{job_id}"
            
            max_attempts = 240
            for attempt in range(max_attempts):
                res = requests.get(status_url)
                status = res.text.strip()
                
                if status == "FINISHED":
                    break
                elif status in ["RUNNING", "PENDING", "QUEUED", "NOT_FOUND"]:
                    time.sleep(5)
                    sys.stdout.write(".")
                    sys.stdout.flush()
                else:
                    print(f"\nError: Job status returned '{status}'")
                    return [], []
            
            if status != "FINISHED":
                print(f"\nError: Search timed out after {max_attempts * 5 / 60} minutes. Status is still: {status}")
                return [], []
                    
            print("\nJob completed. Fetching standard text results...")
            
            # Fetch the standard plaintext output
            out_url = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{job_id}/out"
            out_res = requests.get(out_url)
            out_res.raise_for_status()
            
            homolog_ids = []
            scores = []
            
            lines = out_res.text.split('\n')
            parsing_hits = False
            
            for line in lines:
                # Trigger when we reach the summary table
                if "Sequences producing significant alignments:" in line:
                    parsing_hits = True
                    continue
                
                if parsing_hits:
                    stripped = line.strip()
                    if not stripped:
                        # If we hit a blank line after collecting hits, the table is over
                        if len(homolog_ids) > 0:
                            break
                        else:
                            continue
                            
                    # Example line: "SP:A0A0K8P6T7 PETH_PISS1 Poly(ethylene ...  581     0.0"
                    parts = stripped.split()
                    if len(parts) >= 3:
                        raw_id = parts[0]
                        
                        # Isolate the core UniProt accession
                        if ":" in raw_id:
                            acc = raw_id.split(":")[-1]
                        elif "|" in raw_id:
                            acc = raw_id.split("|")[1]
                        else:
                            acc = raw_id
                        
                        acc = acc.split(".")[0]
                        
                        try:
                            score = float(parts[-2])
                        except ValueError:
                            score = 0.0
                            
                        homolog_ids.append(acc)
                        scores.append(score)
                        
                        if len(homolog_ids) >= hit_limit:
                            break
                            
            print(f"Discovered {len(homolog_ids)} UniProt homologs.")
            return homolog_ids, scores
            
        except Exception as e:
            print(f"\nSearch failed: {e}")
            return [], []