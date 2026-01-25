from Bio.Blast import NCBIWWW, NCBIXML

class DiscoveryEngine:
    def find_homologs(self, sequence_fasta, hit_limit=50):
        print(f"Sending sequence to NCBI BLASTp (be patient)...")
        try:
            # Main BLAST search for initial homologs, lots of room here for optimization/tweaking
            result_handle = NCBIWWW.qblast("blastp", 
                                           "nr", 
                                           sequence_fasta, 
                                           hitlist_size=hit_limit)
            
            blast_record = NCBIXML.read(result_handle)
            homolog_ids = []
            
            for alignment in blast_record.alignments:
                # We take the accession (e.g., 'WP_012345.1')
                accession = alignment.accession
                homolog_ids.append(accession)
                
            print(f"Found {len(homolog_ids)} homologs in the diversity space.")
            return homolog_ids
        except Exception as e:
            print(f" BLAST failed: {e}")
            return []