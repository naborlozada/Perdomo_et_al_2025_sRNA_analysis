#!/bin/bash

# author: Alejandro Nabor Lozada Chavez
# questions: nabor.lozada@gmail.com

# --------------------------------------------------------------------------------------------------------
# Description
# This script process the raw reads of a single sample (R1) by: mapping to a reference sequence, sort it, and calculate (R script call):
#    coverage, size distribution, and small RNAs presence across the targeted reference sequence.
#
# Run: 
#      bash  main_pipeline.sh  [ARG_1]  [ARG_2]
# where: 
#     ARG_1: sample names with the full path directory
#     ARG_2: targeted sequence ID: EVE or CFAV
# --------------------------------------------------------------------------------------------------------


# timing
START_TIME=$(date +%s)


# set workdir
workdir=/home/nlozada/perdomo_project/scripts/
cd workdir;

# arguments:
SAMPLE=$1;
TARGET_REF_SEQ=$2;
FILENAME=$(basename ${SAMPLE} .fq.gz);
RSCRIPT_COVERAGE=sRNAs_analysis_coverage_and_sizes.R;
OUTDIR="/home/nlozada/perdomo_project/outputs/mapping_to_CFAV";
REFERENCE_SEQ2MAP="/home/nlozada/perdomo_project/reference_sequences/CFAV_genome/CFAV_genome_RioPiedras";       ## EVE sequence or CFAV genome



echo
echo "Arguments..."
echo
echo "sample:      $SAMPLE";
echo "filename:    $FILENAME";
echo "RefSeqID:    $TARGET_REF_SEQ";
echo "outdir:      $OUTDIR";
echo "RefSeq2Map:  $REFERENCE_SEQ2MAP";
echo "Rscript:     $make_coverage_analysis_EVE_R";
echo
echo

wait;
sleep 2;



#### ** Main Pipeline Script ** ####

# Alignment
echo "----------------------------------------------------------------------------------------------------------------------------------"
echo "[1] ---> Alignment_command_line: time bowtie -q -v 1 -S -p 20 -x $REFERENCE_SEQ2MAP $SAMPLE  $OUTDIR/${FILENAME}.sam"
# command line
time bowtie -q -v 1 -S -p 20 -x $REFERENCE_SEQ2MAP $SAMPLE  $OUTDIR/${FILENAME}.sam
wait;
sleep 3;


# Sort alignment
echo "----------------------------------------------------------------------------------------------------------------------------------"
echo "[2] ---> Sort_ALN_command_line: time samtools sort  --threads 10  $OUTDIR/${FILENAME}.sam  -o $OUTDIR/${FILENAME}.sorted.sam"
# command line
time samtools sort  --threads 10  $OUTDIR/${FILENAME}.sam  -o $OUTDIR/${FILENAME}.sorted.sam

wait;
sleep 3;


## Calculate coverage, size groups and distribution across the target mapped sequence (R script)
echo "----------------------------------------------------------------------------------------------------------------------------------"
echo "[3] ---> Rscript_command_line: time /home/nlozada/R/R-3.6.2/lib/R/bin/Rscript  $RSCRIPT_COVERAGE  $OUTDIR/${FILENAME}.sorted.sam  ${FILENAME}  $TARGET_REF_SEQ $OUTDIR";
# command line
time /home/nlozada/R/R-3.6.2/lib/R/bin/Rscript  $make_coverage_analysis_EVE_R  $OUTDIR/${FILENAME}.sorted.sam  ${FILENAME}  $TARGET_REF_SEQ $OUTDIR
wait;
sleep 3;

echo
echo "----------------------------------------------------------------------------------------------------------------------------------"
echo


END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

echo "Job complete."
echo "Start Time: $(date -d @$START_TIME)"
echo "End Time: $(date -d @$END_TIME)"
echo "Total Time: ${TOTAL_TIME} seconds"
echo
echo
