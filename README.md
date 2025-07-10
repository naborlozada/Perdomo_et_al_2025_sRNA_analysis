# Perdomo_et_al_2025__Aalbo_CFAV_smallRNAs


Research paper: **Metabolic homeostasis favors tolerance to persistent Cell-fusing agent virus infection in *Aedes aegypti* mosquitoes*** \*

Authors: Hugo Perdomo, Ayda Khorramnejad, Francesco Chen, Alejandro Nabor Lozada-Chávez, Mariangela Bonizzoni.

\* This manuscript is under review.

---

Description:\
This repository contains the datasets (link), and bioinformatics protocol steps used to analyzed the RNAseq data to detect small RNAs in *Aedes albopictus* mosquitoes infected with Cell Fusion Agents (CFAV). 

- datasets:
    **RNAseq, raw reads:** [ONLINE] Original raw reads used for this research paper are publicly available (after publication) at the SRA database under the BioProject ID: **PRJNA1287247**.
       
- bioinformatics pipeline:
   * **BASH command line scripts** used to process the raw reads for quality control and mapping to target sequences (protocol file in Markdown format).
   * **R script** to process and clean mapped reads, find small RNAs, and cuantify their  sequence sizes, coverage, and their distribution across targeted sequences.
   * **MultiQC quality report**. An extended report of the quality control of reads and alignments of each sample.

** Each script has notes with general and/or technical information. For details, read the Methods section in the manuscript. 

