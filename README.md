# Perdomo_et_al_2025__Aalbo_CFAV_smallRNAs


Research paper: **TITLE_HERE*** \*

Running title: "Aedes albopictus, CFAV and small RNAs"

Authors: Hugo Perdomo, Alejandro Nabor Lozada-Chávez, ... , Mariangela Bonizzoni.

\* This manuscript is under review.

---

Description:\
This repository contains the datasets (link), and bioinformatics protocol steps used to analyzed the RNAseq data to detect small RNAs in *Aedes albopictus* mosquitoes infected with Cell Fusion Agents (CFAV). 

- datasets:
    **RNAseq, raw reads:** Original raw reads used for this research paper are publicly available at the SRA database under the ID: *SRA_XXXXX*.
       
- bioinformatics protocol steps:
    1. **BASH command line scripts** used to process the raw reads for quality control and mapping to target sequences (protocol file in Markdown format).
    2. **R script** to process and clean mapped reads, find small RNAs, and cuantify their  sequence sizes, coverage, and their distribution across targeted sequences.
    3. **MultiQC quality report**. An extended report of the quality of reads and alignments of each sample.

** Each script has reference notes with general and/or technical information. For details, read the Methods section in the manuscript. 

