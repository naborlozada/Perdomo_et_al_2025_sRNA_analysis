## Description of the pipeline.

The pipeline was applied to all mosquito samples described in the paper. Raw reads are publicly available in the NCBI SRA database under the BioProject ID PRJNA1287247.
Reads were mapped against two different reference sequences: 1) [Cell Fusion Agent Virus (CFAV) genome "Rio_Piedras GQ165810"](https://www.ncbi.nlm.nih.gov/nuccore/GQ165810.1), 2) [CFAV-EVE-3 from Crava et al., 2020](https://onlinelibrary.wiley.com/doi/10.1111/mec.15798), and 3) EVE-2 from Suziki et al. (2020). See methods in the paper.

A short simple pipeline was created to combine all needed tasks: identify, measure sizes, quantify and locate the distribution of small RNAs of different lengths in the reference sequences separately.

To run properly the pipeline `main_pipeline.sh`, script path files and programs/packages must be previously installed. Additional modifications in the R script should be included in the final output directory (`make perdomo_etal_2025`).

PIPELINE:

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
# CFAV    = Cell Fusion Agent Virus (CFAV) genome "Rio_Piedras GQ165810
# EVE     =  CFAV-EVE-3
# EVE2    =  EVE-2

# built indexes:
bowtie-build  $CFAV  CFAVgenome
bowtie-build  $EVE   EVE2seq
bowtie-build  $EVE2  EVE3seq

# After this, modify the "main_pipeline.sh" script to add each reference sequence, use the index_name (i.e. CFAVgenome)
```

3. Run pipeline:
```bash
# single sample jobs run
for sample in /path/to/raw_reads/mapping_to_CFAV/*.fq.gz; do bash main_pipeline.sh  $sample CFAV  2>&1 | tee ${sample}.CFAV.stderr.log; 
for sample in /path/to/raw_reads/mapping_to_EVE/*.fq.gz; do bash main_pipeline.sh  $sample EVE  2>&1 | tee ${sample}.EVE.stderr.log;
```



