import Bio
from Bio import Entrez
import pandas as pd

print (f" Environment Active: Biopython {Bio.__version__}")
print (f" Pandas {pd.__version__} is ready for dataframes.")

# Check if we can talk to NCBI
Entrez.email = "aangushenry@gmail.com"
print ("Testing connection to NCBI...")

handle = Entrez.esearch(db="nucleotide", term="Methyl parathion hydrolase", retmax=1)

record = Entrez.read(handle)

print(f"Connection successful! Found {record['Count']} matches.")

