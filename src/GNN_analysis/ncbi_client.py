import os
import time
import socket
from Bio import Entrez, SeqIO
from dotenv import load_dotenv

# Set global socket timeout to prevent stalled data from ncbi (occurs frequently)
socket.setdefaulttimeout(65)

load_dotenv()
Entrez.email = os.getenv("NCBI_EMAIL")
Entrez.api_key = os.getenv("NCBI_API_KEY")

class NCBIClient:
    def get_locations(self, protein_id):
        # First loop: IPG database
        for attempt in range(4):
            try:
                locations = self._get_ipg_locations(protein_id)
                if locations:
                    return locations
                
                else:
                    # Server successful but no IPG record
                    # Switch immedately to GenPept
                    print(f"    -> No IPG data found for {protein_id}, switching to GenPept")
                    break

            # Exception for first loop IPG 
            except Exception as e:
                if attempt < 3: 
                    wait_time = 2 ** attempt
                    print(f"    -> Failed to fetch IPG report for {protein_id}, retrying in {wait_time}s... ")
                    time.sleep(wait_time)
                else:
                    print(f"    -> Failed to fetch IPG report for {protein_id}, switching to GenPept")

        # Second loop: GenPept database
        for attempt in range (4):
            try: 
                locations = self._get_direct_locations(protein_id)
                if locations:
                    return locations
                
                else:
                    # Server successful but no GenPept report
                    break

            # Second loop exception    
            except Exception as e:
                if attempt < 3:
                    wait_time = 2 ** attempt
                    print(f"    -> Failed to fetch GenPept data for {protein_id}, retrying in {wait_time}s...")
                    time.sleep(wait_time)
                else:
                    print(f"    -> Failed to fetch GenPept data for {protein_id}, skipping...")

        # Fallback for total failure of protein
        return[]

    def _get_ipg_locations(self, protein_id):
            handle = None
            try:
                handle = Entrez.efetch(db="protein", id=protein_id, rettype="ipg", retmode="xml", timeout=45)
                raw_data = Entrez.read(handle)
                locations = []
                if isinstance(raw_data, dict) and 'IPGReport' in raw_data:
                    report = raw_data['IPGReport']
                    if isinstance(report, list):
                        for protein in report: self._extract_cds_ipg(protein, locations)
                    else:
                        self._extract_cds_ipg(report, locations)
                return locations
            finally: 
                if handle:
                    handle.close()

    def _get_direct_locations(self, protein_id):
        handle = None
        try:
            handle = Entrez.efetch(db="protein", id=protein_id, rettype="gp", retmode="text", timeout=45)
            record = SeqIO.read(handle, "genbank")
            locations = []
            for feat in record.features:
                if feat.type == "CDS":
                    coded_by = feat.qualifiers.get("coded_by", [""])[0]
                    if coded_by:
                        parts = coded_by.split(":")
                        nuc_acc = parts[0]
                        strand = -1 if "complement" in coded_by else 1
                        coord_str = parts[1].replace("complement(", "").replace(")", "").replace("<", "").replace(">", "")
                        coords = coord_str.split("..")
                        locations.append({'nuc_acc': nuc_acc, 'start': int(coords[0]), 'end': int(coords[1]), 'strand': strand})
            return locations
        finally: 
            if handle:
                handle.close()

    def _extract_cds_ipg(self, protein_data, locations_list):
        if not isinstance(protein_data, dict): return
        cds_list = protein_data.get('CDSList', [])
        for cds in cds_list:
            attr = cds.attributes if hasattr(cds, 'attributes') else cds
            locations_list.append({
                'nuc_acc': attr.get('accver'), 
                'start': int(attr.get('start', 0)), 
                'end': int(attr.get('stop', 0)), 
                'strand': 1 if attr.get('strand') == '+' else -1
            })

    def fetch_neighborhood(self, nuc_acc, start, end, window=10000):
        f_start, f_end = max(1, start - window), end + window
        # Added a small retry loop for network stability
        
        
        for attempt in range(4):
            handle = None
            try:
                # API rate limit buffer
                time.sleep(0.15)

                handle = Entrez.efetch(db="nuccore", id=nuc_acc, seq_start=f_start, seq_stop=f_end, rettype="gbwithparts", retmode="text", timeout=60)
                record = SeqIO.read(handle, "genbank")
                return record
            
            except Exception as e:
                if attempt < 3 :
                    wait_time = 2 ** attempt
                    print(f"Network error on {nuc_acc}, retrying in {wait_time}s...")
                    time.sleep(wait_time)
                else:
                    print(f"Timeout fetching neighborhood for {nuc_acc} after 4 attempts. Skipping.")
            finally:
                if handle:
                    handle.close()
        return None