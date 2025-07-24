## Description of the pipeline.

The pipeline was applied to all mosquito samples described in the paper. Raw reads are publicly available in the NCBI SRA database under the BioProject ID PRJNA1287247.
Reads were mapped against three different reference sequences independently (for details see methods in the paper): 
  1) [Cell Fusion Agent Virus (CFAV) genome "Rio_Piedras GQ165810"](https://www.ncbi.nlm.nih.gov/nuccore/GQ165810.1)
  2) [CFAV-EVE-3 from Crava et al., 2020](https://onlinelibrary.wiley.com/doi/10.1111/mec.15798)
  3) EVE-2 from Suziki et al. (2020). 

A short simple pipeline was created to combine all needed tasks: identify, measure sizes, quantify and locate the distribution of small RNAs of different lengths in the reference sequences separately. See notes below.

NOTE 1: For simplicity, the pipeline `main_pipeline.sh` script is writen to run only one reference sequence at time, that is, only results associated to the reference sequences `CFAV_genome_RioPiedras` will be produced (after adapting the line of code to your local path and retaining the same index name). To get the results for the other two reference sequences (EVE2 and EVE3) replace the index name and full path file with the corresponding reference sequence. Line to modify: `REFERENCE_SEQ2MAP`.

NOTE 2: To run properly the pipeline `main_pipeline.sh` script, full path files (reference sequence infiles, output directories) and programs/packages (binaries, R) must be previously installed and set it properly. Programs used: bowtie, samtools, R, and Tidyverse, Biostrings and ggplot2 R packages. See details in Methods in the manuscript. The location of the R script should be declare in the variable `RSCRIPT_COVERAGE` of the `main_pipeline.sh` script. 

START.

1. Quality Control analysis.
```bash
# make stats
for sample in /path/to/raw_reads/mapping_to_*/*.fq.gz; do samtools flagstats $sample > ${sample}.txt;
for sample in /path/to/raw_reads/mapping_to_*/*.fq.gz; do fastqc -o outputdir/ $sample;

# summarize stats:
cd mapping_to_CFAV/
multiqc .
cd ../mapping_to_EVE/
multiqc .
 ```


2. Index target reference sequences:
```bash
# Reference sequences
CFAV=/path/to/refseq/CFAV/CFAV_genome_Rio_Piedras    # Cell Fusion Agent Virus (CFAV) genome "Rio_Piedras GQ165810
EVE3=/path/to/refseq/EVE3/EVE3                       # CFAV-EVE-3
EVE2=/path/to/refseq/EVE2/EVE2                       # EVE-2

# syntaxis:
#  program  fasta  index_name

# built indexes:
bowtie-build  ${CFAV}.fasta  CFAVgenome
bowtie-build  ${EVE}.fasta   EVE2seq
bowtie-build  ${EVE2}.fasta  EVE3seq

# After this, modify the "main_pipeline.sh" script to add each reference sequence (as described aboved). Be sure to use the index_name (i.e. CFAVgenome)
```


3. Run pipeline:
```bash
# Single sample job run (assuming there is not much computational resources). To reduce/increase the computation power,
# change the value in the options "-p" of bowtie (now set it to use 20 threads), or the option "thread" in samtools. The R script runs in a single core.
# Overall, it uses around 30Gb of RAM.
 
for sample in /path/to/raw_reads/mapping_to_CFAV/*.fq.gz; do bash main_pipeline.sh  $sample CFAV  2>&1 | tee ${sample}.CFAV.stderr.log; 
for sample in /path/to/raw_reads/mapping_to_EVE/*.fq.gz; do bash main_pipeline.sh  $sample EVE  2>&1 | tee ${sample}.EVE.stderr.log;
```



