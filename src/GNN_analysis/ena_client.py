import aiohttp
import asyncio
from Bio import SeqIO
import io

class UniProtENAClient:
    def __init__(self, max_concurrent=15):
        # Prevent hammering the EBI/UniProt APIs
        self.semaphore = asyncio.Semaphore(max_concurrent)

    async def fetch_uniprot_data(self, session, uniprot_id):
        """Fetches both mapping locations and rich pathway/domain metadata."""
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        
        async with self.semaphore:
            try:
                async with session.get(url, timeout=20) as response:
                    if response.status != 200:
                        return None, []
                    
                    data = await response.json()
                    
                    metadata = {
                        "review_status": data.get("entryType", "Unknown"),
                        "ec_number": self._extract_ec(data),
                        "interpro_ids": self._extract_xref(data, "InterPro"),
                        "pfam_ids": self._extract_xref(data, "Pfam"),
                        "go_terms": self._extract_xref(data, "GO"),
                        "pathway_xrefs": self._extract_xref(data, "PathwayCommons"),
                        "has_alphafold": 1 if self._extract_xref(data, "AlphaFoldDB") else 0,
                        "sequence": data.get("sequence", {}).get("value", "")
                    }
                    
                    locations = self._extract_embl_locations(data)
                    return metadata, locations
            except (asyncio.TimeoutError, aiohttp.ClientError):
                return None, []

    def _extract_ec(self, data):
        try:
            return data['proteinDescription']['recommendedName']['ecNumbers'][0]['value']
        except KeyError:
            return ""

    def _extract_xref(self, data, db_name):
        xrefs = data.get("uniProtKBCrossReferences", [])
        ids = [ref["id"] for ref in xrefs if ref["database"] == db_name]
        return ", ".join(ids)

    def _extract_embl_locations(self, data):
        locations = []
        xrefs = data.get("uniProtKBCrossReferences", [])
        for ref in xrefs:
            if ref["database"] == "EMBL":
                embl_acc = ref["id"]
                locations.append({"embl_acc": embl_acc})
        return locations

    async def fetch_ena_neighborhood(self, session, embl_acc):
        """Fetches the EMBL flat file from ENA asynchronously."""
        url = f"https://www.ebi.ac.uk/ena/browser/api/embl/{embl_acc}?download=false"
        
        async with self.semaphore:
            try:
                # ENA API can be slow; allow slightly longer timeout
                async with session.get(url, timeout=45) as response:
                    if response.status != 200:
                        return None
                    text_data = await response.text()
                    
                    try:
                        # Catch Exception rather than ValueError to absorb Biopython AssertionErrors
                        record = SeqIO.read(io.StringIO(text_data), "embl")
                        return record
                    except Exception:
                        return None
                        
            except (asyncio.TimeoutError, aiohttp.ClientError):
                return None