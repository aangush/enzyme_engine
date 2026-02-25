# Enzyme Engine

### A Comparative Genomics Pipeline for Discovering Functional Gene Modules and Exploring Bacterial Protein Diversity

#### **Author: Aidan William Angus-Henry**

---

## Overview

Enzyme Engine is a bioinformatics tool designed to find and map genomic neighborhood networks of prokaryotic protein homologs. It combines BLAST discovery with local synteny analysis and domain identification to reveal conserved genetic modules across phyla. Features in-development include known operon function mapping through integration with MIBiG and MetaCyc databases, followed by a "deep search" consisting of phylogenetic and structural analysis of un-annotated/hypothetial proteins of interest.

### Motivation

This project was inspired by the technical limitations and bottlenecks I encountered working with poorly characterized biochemical pathways for metabolic engineering, and by my use of the [Enzyme Function Initiative Tools](https://efi.igb.illinois.edu/). It is an attempt to streamline, optimize, and expand on comparative genomic context searching to discover functional "accessory" proteins for a pathway of interest.

_Furthermore, this project was motivated by a desire to leverage comparative genomic and environmental context to characterize and harness biodiversity, offering a targeted alternative to black-box predictive models trained on insufficient data._

---

## Key Features

- Takes a FASTA formatted protein sequence, and runs a large BLASTp optimized to find homologs (="anchors") and perform a representative search through bacterial biodiversity.
- Uses signed genomic distances to distinguish upstream and downstream neighbors of anchor proteins, and reveal possible operon orientation.
- Automatically runs `hmmscan` against a Pfam-A database (locally) to annotate domains in hypothetical proteins.
- Calculates and displays frequency, avg distance, spread of avg distances, strand, and avg size of genomic neighbors of homologs (="anchor" proteins).
- Finds and reports top co-occurring genomic neighbors and modules for follow-up analysis.

---

## Getting Started

### Prerequisites

- Python 3.10
- HMMER3 Installed and accessible in PATH
- NCBI API Key (Highly recommended for 100+ hits)
- Local Pfam-A Database formatted with `hmmpress`

### Installation
1. Clone the repository and enter the directory:
`git clone [https://github.com/aangush/enzyme_engine.git](https://github.com/aangush/enzyme_engine.git)
cd enzyme_engine`

2. Create and activate the Conda environment:
`conda env create -f environment.yml 
conda activate enzyme_env`

3. Set up environmental variables: Create a .env file in the root directory and add your NCBI credentials to prevent strict API limiting
`NCBI_EMAIL=your.email@domain.com
NCBI_API_KEY=your_api_key_here`

4. Prepare the pfam database: Download the `Pfam-A.hmm` database, format it using `hmmpress Pfam-A.hmm`, and place all resulting files into a data/pfam directory within the project root.

### Quick Start:
To run the pipeline, provide a FASTA file containing your target anchor protein sequence. Execute the main script from the root directory:
`python src/GNN_analysis/main.py path/to/input.fasta`

---


### Workflow

1. Performs a BLASTp search to identify N homologs of input sequence.
2. Queries NCBI for the genomic location of each homolog, and fetches a +-10kb window.
3. Scans every hypothetical neighbor against Pfam and cluster remaining hypothetical proteins by sequence identity.
4. Reports a synteny summary with frequency, avg distance, spread of distances, strand, and avg size of neighbors for easy viewing.
5. Reports a top co-occurring neighbor pair table for given neighbors A and B with frequencies of co-occurrence.
6. Runs and reports a co-occurrence analysis searching for linked modules of multiple members using genomic neighborhood data.

---

### Example output with input PETase (A0A0K8P6T7.1) plastic degrading enzyme from _Piscinibacter sakaiensis_ for 200 BLAST hits:

![example enzyme engine output using PETase protein sequence as input](docs/PETase_output_1.png)

**Genomic Neighbors Overview**

- Top genomic neighbors of PETase homologs include chaperones, transporters, other hydrolases, and efflux pumps.
- Spread is a function of the (maximum largest distance away - minimum largest distance away), providing a crude estimate of variation of nieghbor genomic distance to anchor.

![example enzyme engine functional module suggestion output with PETase input](docs/PETase_output_modules.png)

**Co-occurrence analysis**

- Enzyme Engine displays the top 10 functional modules found sorted by frequency with at least 4 members, providing insights about functional linkage.
- For our 200-hit PETase search, the most frequent module contains a full transporter, efflux pump, and dihydroorotase. Whether or not dihydroorotase has a function in microbial PET degredation remains unknown.

---

### Known Limitations & Technical Notes
**NCBI Data Attrition and API Bottlenecks:** Due to the highly fragmented and uncurated nature of environmental sequence data (e.g., short WGS contigs missing flanking regions), NCBI Entrez queries for genomic neighborhoods frequently return empty. Furthermore, the NCBI API is prone to transient stalls and timeouts when fetching bulk GenBank records. As a result, a percentage of initial BLAST hits will inevitably be dropped during the neighborhood retrieval phase. I am currently working on scaling the query architecture  to process much larger initial hit pools, which will compensate for this data attrition and ensure high statistical power for functional module detection (and also take much longer).

---


## Future Directions and Features (In development)

- Increasing BLASTp hit number and breadth to better probe genomic space and gather more data.
- Expanding multiple-linkage analysis to probe larger genomic windows than the current +-10kb.
- Pull existing metabolic pathway data from MetaCyc and assign putative operon members to precise biochemical reactions to find mystery proteins that seem to be linked with a given gene cluster.
- For a given hypothetical/un-annotated protein, trigger a "deep search" that integrates phylogenetic, biosample, and structural analysis to investigate functional role of linked protein.
- Alphafold integration, active site analysis, and possible molecular dynamics simulations to probe for functional mutations or sufficient deviation from predicted function.
