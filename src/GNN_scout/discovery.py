from Bio.Blast import NCBIWWW, NCBIXML

class DiscoveryEngine:
    def find_homologs(self, sequence_fasta, hit_limit=200):
        print(f"Sending sequence to NCBI BLASTp (be patient)...")
        try:
            # Main BLAST search for initial homologs, lots of room here for optimization/tweaking
            result_handle = NCBIWWW.qblast("blastp", 
                                           "nr", 
                                           sequence_fasta, 
                                           hitlist_size=hit_limit)
            
            blast_record = NCBIXML.read(result_handle)
            homolog_ids = []
            identities = []
            
            for alignment in blast_record.alignments:
                # Take the accession (e.g., 'WP_012345.1')
                accession = alignment.accession
                homolog_ids.append(accession)

                # Get the first HSP for alignment
                if alignment.hsps:
                    hsp = alignment.hsps[0]

                    # Identities and alignment length are attributes of HSP
                    ident = (hsp.identities / hsp.align_length) * 100
                    identities.append(ident)
                else:
                    identities.append(0)
                
            print(f"Found {len(homolog_ids)} homologs in the diversity space.")
            
            
            return homolog_ids, identities
        
        except Exception as e:
            print(f" BLAST failed: {e}")
            return [], []