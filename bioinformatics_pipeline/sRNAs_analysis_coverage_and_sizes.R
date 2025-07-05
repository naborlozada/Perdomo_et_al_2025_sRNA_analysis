#!/usr/bin/env/R



# Author: Alejandro Nabor Lozada Chávez
# Questions: nabor.lozada@gmail.com


# ---------------------------------------------------------------------------------------------------------------------------------------
# Description:
# ---------------------------------------------------------------------------------------------------------------------------------------
# R scripts that quantifies and obtains the distribution of sRNAs when targeted to a specific reference genome or sequence.
#
# The distribution and quantification of small RNAs across for a reference sequence it is estimated with this custom R script using
# simple code and functions to parse a SAM file with the tidyverse v1.3.1 and Biostrings v2.54.0 R packages.
#
# Briefly, sequence length, mapping strands and reads counts information were extracted from SAM files, then small RNAs count frequencies
# mapped either to one or both strands were performed only for sequence lengths ≥18 and ≤33 bp. Outputs are produced in CSV format files.
# Figures were create with PRISM GradPath using the CSV files.  
# ---------------------------------------------------------------------------------------------------------------------------------------


# Expected time with an EVE sequence:
#    user  system elapsed 
# 826.032  42.212 873.693 





suppressMessages(library(tidyverse));
suppressMessages(library(ggplot2));
suppressMessages(library(Biostrings));




rm(list=ls());



# Assign arguments to variables
args        <- commandArgs(trailingOnly = TRUE);
sam_infile  <- args[[1]];
sample      <- args[[2]];
treatment   <- args[[3]];
output_path <- args[[4]];


# define length of reference genome:
seqTargetLength = 0;

if (treatment=="EVE") {
    seqTargetLength = 734;
    #print("it is positive.")
} else if (treatment=="CFAV") {
    seqTargetLength = 9765;
} else {
    cat("\n")
    print("WARNING: 'seqTargetLength' NOT ASSIGNED! Please, define the target mapped sequence (3rd arg): 'CFAV' OR 'EVE'.")
    stop("The value is ZERO, so the script must end here.")
}



cat("\n\nInput arguments:\n")
sam_infile
sample
treatment
output_path
seqTargetLength
cat("\n\n")




# set dir
setwd(output_path);


# set main filename
sample_name <-  paste0(sample,".",treatment);



# load full SAM FILE
cat("loading SAM file (time)...\n");
# ---------------------------------------------------------------------------------------------------------------------------------------
system.time( SAM_file <- read.delim(sam_infile, header = FALSE,row.names = NULL, sep = "")[,1:12] );

cat("\n\nContinue with parsing data...")


# parse data
SAM_file_cleaned_A <- SAM_file %>% 
                        dplyr::mutate_if(is.factor, as.character) %>% 
                        dplyr::glimpse()



# before anything else, remove SAM headers
SAM_file_cleaned_B <- SAM_file_cleaned_A %>% dplyr::filter(!grepl("^@", V1))

# header names
colnames(SAM_file_cleaned_B) <- c("QNAME", "FLAG", "RNAME", "POS", "MAPQ","CIGAR", "MRNM", "MPOS","ISIZE", "SEQ", "QUAL", "OPT");


# change col properties
SAM_file_cleaned_B$FLAG   <- as.numeric(as.character(SAM_file_cleaned_B$FLAG));
SAM_file_cleaned_B$POS    <- as.numeric(as.character(SAM_file_cleaned_B$POS));
SAM_file_cleaned_B$MAPQ   <- as.numeric(as.character(SAM_file_cleaned_B$MAPQ));
SAM_file_cleaned_B$MPOS   <- as.numeric(as.character(SAM_file_cleaned_B$MPOS));
SAM_file_cleaned_B$ISIZE  <- as.numeric(as.character(SAM_file_cleaned_B$ISIZE));



# add some new information
sample_data_reads <- mutate(SAM_file_cleaned_B,  SAMPLE = sample, Treatment = treatment, SEQ_RC = chartr("ATGC", "TACG", SEQ), alignmentLength = nchar(as.vector(SEQ)));





cat("Make counts for LENGTH SEQUENCE size distribution...")
# ---------------------------------------------------------------------------------------------------------------------------------------
# total number of reads:
total_reads_count <- nrow(sample_data_reads)

WARNING_MSM = "WARNING";


# Parse sense/antisense strand and make counts based on length reads 
sample_data_reads_stats <- dplyr::filter(sample_data_reads) %>%
                                dplyr::select(Treatment, alignmentLength, FLAG) %>%
                                        dplyr::group_by(Treatment, FLAG, alignmentLength) %>%
                                            dplyr::summarise(Counts = n()) %>%
                                                dplyr::mutate(ReadsCounts_byStrand = ifelse(FLAG == 16, "-",
                                                                                     ifelse(FLAG == 0, "+", "unmapped"))
                                                             ) %>% as.data.frame();



# remove unmapped (if there is any)
sample_data_reads_stats.flitrd <- dplyr::filter(sample_data_reads_stats, FLAG != 4) %>% 
                                                dplyr::select(alignmentLength, ReadsCounts_byStrand, Counts)



# Create empty sized dataframe

# variables
qwidth_values <- rep(18:33, each = 2)
strand_values <- rep(c("+", "-"), times = 16)
frequency_values <- rep(0, 16)

# table
data2fill <- data.frame(alignmentLength = qwidth_values, 
                 ReadsCounts_byStrand = as.factor(strand_values), 
                 Counts = frequency_values)


# Merge the new dataframe into the initial dataframe, updating only the matching qwidth values
data2plot <- merge(data2fill, sample_data_reads_stats.flitrd, by = c("alignmentLength", "ReadsCounts_byStrand"), all.x = TRUE)


# Parse the Frequency data values and replace them with the new datapoints
data2plot$Frequency <- ifelse(!is.na(data2plot$Counts.y), 
                                     data2plot$Counts.y, 
                                     data2plot$Counts.x)

# Select columns
data2plot_updated <- data2plot[, c("alignmentLength", "ReadsCounts_byStrand", "Frequency")]



# /// output data table for plot ///
filename_ReadsSizeDistr_CSV <- paste0(output_path,"/",sample_name,".size_distribution.both_strands.csv");

write.csv(data2plot_updated, file=filename_ReadsSizeDistr_CSV, row.names = FALSE)




# make plot: stacked-bar
sample_mapped_reads_both_strands <-  ggplot(data2plot_updated, aes(fill=as.factor(ReadsCounts_byStrand), y=Frequency, x=alignmentLength)) + 
                                        geom_bar(position="dodge", stat="identity") +
                                        theme_classic() +
                                        scale_fill_manual(values = c("+" = "#f5ad3b", "-" = "blue")) +
                                        scale_x_continuous(expand = c(0, 0.1), breaks = seq(18, 33, by = 1)) +  ### limits = c(18,34), 
                                        scale_y_continuous(expand = c(0.005, 0)) +
                                        guides(fill=guide_legend(title="Strand")) +
                                        labs(x = "Size (bp)", y = "Number of reads")


filename_ReadsSizeDistr_PLOT <- paste0(output_path,"/",sample_name,".size_distribution.both_strands.pdf");

pdf(file=filename_ReadsSizeDistr_PLOT, width=15, height=10)
sample_mapped_reads_both_strands
dev.off()



cat("Table and plots are done...\n\n")




cat("Make GENOME COVERAGE data...\n\n")
# ---------------------------------------------------------------------------------------------------------------------------------------
# remove not mapped reads...
sample_data_reads_flitrd <- dplyr::filter(sample_data_reads, FLAG != 4);


# select columns
sample_data_reads_flitrd_simplified <- sample_data_reads_flitrd %>% dplyr::select(RNAME,FLAG,POS,alignmentLength,SAMPLE)



# make counts by groups and define mapped stranded or unmapped reads
sample_data_reads_flitrd_simplified.plot <- sample_data_reads_flitrd_simplified %>%
                                                        dplyr::group_by(FLAG, POS, alignmentLength) %>%
                                                            dplyr::summarise(Count = n()) %>%
                                                                dplyr::mutate(Direction = ifelse(FLAG == 16, Count*-1, Count)) %>%
                                                                dplyr::mutate(ReadsCounts_byStrand = ifelse(FLAG == 16, "-",
                                                                                                     ifelse(FLAG == 0, "+", "unmapped"))
                                                                )  %>% as.data.frame()



# filter data to plot
sample_data_reads_fltrd2plot <- sample_data_reads_flitrd_simplified.plot %>% dplyr::select(POS, ReadsCounts_byStrand, Count, alignmentLength) 

# change data col format 
sample_data_reads_fltrd2plot$ReadsCounts_byStrand <- as.factor(as.character(sample_data_reads_fltrd2plot$ReadsCounts_byStrand));




# -----------------------------------------------------------------------------
# CREATE MAIN TABLE WITH ALL GENOME POSITIONS
# -----------------------------------------------------------------------------

# Create empty dataframe
POS                  <- rep(1:seqTargetLength, each = 2);
ReadsCounts_byStrand <- rep(c("+", "-"), times = seqTargetLength);
Count                <- rep(0, seqTargetLength);
alignmentLength      <- rep(0, seqTargetLength);


cov_data2fill <- data.frame(POS = POS, 
                            ReadsCounts_byStrand = as.factor(ReadsCounts_byStrand), 
                            Count = Count,
                            alignmentLength = alignmentLength)


# merge both data tables
# -----------------------------------------------------------------------------
# Merge the new dataframe into the initial dataframe, updating only the matching qwidth values
coverage_data2plot <- base::merge(cov_data2fill, sample_data_reads_fltrd2plot, by = c("POS", "ReadsCounts_byStrand", "alignmentLength"), all = TRUE)




# Update the Frequency values in the initial dataframe with those from the new dataframe
coverage_data2plot$Frequency <- ifelse(!is.na(coverage_data2plot$Count.y), 
                                              coverage_data2plot$Count.y, 
                                              coverage_data2plot$Count.x)



cat("\n\n\n")
cat("COVERAGE DATA 'coverage_data2plot': top 100 lines")
cat("\n------------------------------------------------------------------------------------------------------------------------\n")
head(coverage_data2plot,n=100)
cat("\n------------------------------------------------------------------------------------------------------------------------\n")
cat("\n\n")

# Select relevant columns and rename them to the original column names
coverage_data2plot_updated <- coverage_data2plot[, c("POS","alignmentLength", "ReadsCounts_byStrand", "Frequency")]


coverage_data2plot_updated <- coverage_data2plot_updated %>% 
                                    dplyr::mutate(Frequency_byStrand = ifelse(ReadsCounts_byStrand == "-", Frequency*-1, Frequency)) %>%
                                    as.data.frame()


# Rename cols
colnames(coverage_data2plot_updated) <- c("Postition","alignmentLength", "Strand", "Frequency", "Frequency_byStrand");
# -----------------------------------------------------------------------------





# WHOLE EVE INTEGRATION SEQUENCE COVERAGE: ALL SEQ-NTS SIZES  
# -----------------------------------------------------------------------------
# output data table for plot
filename_Coverage_allREADS_CSV <- paste0(output_path,"/",sample_name,".coverage.sRNA_sequence_all_reads.csv");

write.csv(coverage_data2plot_updated, file=filename_Coverage_allREADS_CSV, row.names = FALSE)




# NOTE: If required, use this table to classify sRNAs as 21 nts or those associated to a length between 26-to-30 nts.

# save this session
filename_Rsession <- paste0(output_path,"/Perdomo_etal_2025.small_RNAs.seq_size_distribution_and_coverage.",sample_name,".Rsession.RData");

save.image(file = filename_Rsession);



