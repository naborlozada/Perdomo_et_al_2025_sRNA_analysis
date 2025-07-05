## Description of the pipeline.

The pipeline was applied to all mosquito samples described in the paper. Raw reads are publicly available in the NCBI SRA database under the ID SRAXXXXX.
Reads were mapped against two different reference sequences: 1) [Cell Fusion Agent Virus (CFAV) genome "Rio_Piedras GQ165810"](https://www.ncbi.nlm.nih.gov/nuccore/GQ165810.1), and 2) [Endogenous Viral Element 3 from Crava et al., 2020](https://onlinelibrary.wiley.com/doi/10.1111/mec.15798). See methods in the paper.

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

2. Run pipeline:
```bash
for sample in /path/to/raw_reads/mapping_to_CFAV/*.fq.gz; do bash main_pipeline.sh  $sample CFAV  2>&1 | tee ${sample}.CFAV.stderr.log; 
for sample in /path/to/raw_reads/mapping_to_EVE/*.fq.gz; do bash main_pipeline.sh  $sample EVE  2>&1 | tee ${sample}.EVE.stderr.log; 
```



