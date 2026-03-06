import requests
import time
import sys
import io
from Bio.Blast import NCBIXML

class DiscoveryEngine:
    def find_homologs(self, sequence_fasta, hit_limit=500):
        print("Submitting sequence search to EBI Job Dispatcher (NCBI BLAST+ against UniProtKB)...")
        
        # 1. Ensure FASTA formatting
        seq_clean = sequence_fasta.strip()
        if not seq_clean.startswith(">"):
            seq_clean = f">Query_Sequence\n{seq_clean}"

        # 2. EBI Job Dispatcher endpoint and parameters
        submit_url = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run"
        payload = {
            "email": "enzyme_engine_user@gmail.com", # Swapped from @example.com to prevent spam filters
            "program": "blastp",
            "stype": "protein",
            "sequence": seq_clean,
            "database": "uniprotkb"
        }
        
        try:
            # 3. Submit Job
            response = requests.post(submit_url, data=payload)
            if response.status_code != 200:
                print(f"API Error: Server rejected payload. Details: {response.text}")
                return [], []
                
            job_id = response.text.strip()
            print(f"Job submitted successfully. Job ID: {job_id}")
            print("Polling for results (this usually takes 1-3 minutes)...")
            
            status_url = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/{job_id}"
            
            # 4. Poll Status
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
                    
            print("\nJob completed. Fetching results...")
            
            # 5. Try XML first, if 400 Bad Request, fetch the standard output log to diagnose
            result_url = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{job_id}/xml"
            xml_res = requests.get(result_url)
            
            if xml_res.status_code == 400:
                print(f"\nCould not fetch XML (400 Bad Request). Checking standard text output to diagnose...")
                
                # Fetch plain text output
                out_url = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{job_id}/out"
                out_res = requests.get(out_url)
                
                if out_res.status_code == 200:
                    print(f"\n--- SERVER LOG PREVIEW ---\n{out_res.text[:1500]}\n--------------------------")
                    if "No hits found" in out_res.text:
                        print("\nConclusion: The BLAST search completed but found 0 homologs for this sequence.")
                else:
                    # If even the text log is missing, dump the available file types
                    types_url = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/resulttypes/{job_id}"
                    types_res = requests.get(types_url)
                    print(f"Available result files on server: {types_res.text}")
                
                return [], []
                
            xml_res.raise_for_status()
            
            # 6. Parse the XML if successful
            blast_record = NCBIXML.read(io.StringIO(xml_res.text))
            homolog_ids = []
            scores = []
            
            for alignment in blast_record.alignments:
                if len(homolog_ids) >= hit_limit:
                    break
                    
                raw_acc = alignment.accession
                
                if ":" in raw_acc:
                    acc = raw_acc.split(":")[-1]
                elif "|" in raw_acc:
                    acc = raw_acc.split("|")[1]
                else:
                    acc = raw_acc
                    
                homolog_ids.append(acc)
                
                if alignment.hsps:
                    scores.append(alignment.hsps[0].score)
                else:
                    scores.append(0)
                    
            print(f"Discovered {len(homolog_ids)} UniProt homologs.")
            return homolog_ids, scores
            
        except Exception as e:
            print(f"\nSearch failed: {e}")
            return [], []