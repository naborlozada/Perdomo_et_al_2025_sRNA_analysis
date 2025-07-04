#!/usr/bin/env/R



# Author: Alejandro Nabor Lozada Chávez
# Questions: nabor.lozada@gmail.com



# ---------------------------------------------------------------------------------------------------------------------------------------
# Description:
# ---------------------------------------------------------------------------------------------------------------------------------------
#
# R scripts that quantifies and gets the distribution of sequences when targeted to a specific reference genome or sequence.
#
# The distribution and quantification of small RNAs across for a reference sequence it is estimated with this custom R script using
# simple code and functions to parse a SAM file with the tidyverse v1.3.1 and Biostrings v2.54.0 R packages.
#
# Briefly, sequence length, mapping strands and reads counts information were extracted from SAM files, then small RNAs count frequencies
# mapped either to one or both strands were performed only for sequence lengths ≥18 and ≤33 bp.  
#
# ---------------------------------------------------------------------------------------------------------------------------------------


# Expected time with an EVE sequence:
#    user  system elapsed 
# 716.428  40.076 790.099 
#    user  system elapsed 
# 826.032  42.212 873.693 


# VERIFY EVE sequence length:
#   awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' /home/tigerpv_das/users_data/Nabor/otherWorks/hugo/reference_sequences/EVE_integration/EVE_integration.fasta
# >CFAV-EVE-3_viral_integration
# 734
# VERIFY CFAV sequence length:
#   awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' /home/tigerpv_das/users_data/Nabor/otherWorks/hugo/reference_sequences/CFAV_genome/CFAV_genome_RioPiedras.fasta >GQ165810_CFAV_Rio_Piedras
# >GQ165810_CFAV_Rio_Piedras
# 9765
# ---------------------------------------------------------------------------------------------------------------------------------------




suppressMessages(library(tidyverse));
suppressMessages(library(ggplot2));
suppressMessages(library(Biostrings));
#library(data.table)




rm(list=ls());



# Assign arguments to variables
args        <- commandArgs(trailingOnly = TRUE);
sam_infile  <- args[[1]];
sample      <- args[[2]];
treatment   <- args[[3]];
output_path <- args[[4]];


# define length of reference genome:
#ifelse(args[[3]]=="EVE", seqTargetLength = 734, ifelse(args[[3]]=="EVE", seqTargetLength = 9765, "WARNING: 'seqTargetLength' NOT ASSIGNED!!")
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

# sam_infile  =  "L4_1.sorted.sam";
# sample      =  "L4_1";
# treatment   =  "CFAV"
# output_path =  "/home/tigerpv_das/users_data/Nabor/otherWorks/hugo/scripts/main_scripts/final_reproducible_scripts/";
# seqTargetLength =  9765 (CFAV) or 734 (EVE)

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




# stacked-bar plot
sample_cov_plot_both_strands <- ggplot(coverage_data2plot_updated, aes(fill=Strand, y=Frequency_byStrand, x=Postition)) + 
                                    geom_bar(position="dodge", stat="identity", , width=10.0) +
                                    scale_fill_manual(values= c("red", "blue"),
                                                    breaks= c("+", "-"),
                                                    labels= c("+","-")
                                    ) +
                                    scale_x_continuous(expand = c(0, 0.1), limits = c(1,seqTargetLength), breaks = seq(0,seqTargetLength, by = 100)) +
                                    scale_y_continuous(expand = c(0.005, 0), limits = c(-30,80), breaks = seq(-30,80, by = 10)) +
                                    theme_minimal() +
                                    theme(axis.text.x = element_text(angle = 0, hjust = 0.50, vjust = 1.0, colour = "black", size=20),
                                        axis.text.y = element_text(colour = "black",size=20),
                                        axis.title.x = element_text(color="black", size=23, face="bold"),
                                        axis.title.y = element_text(color="black", size=23, face="bold"),
                                        axis.text=element_text(size=10),
                                        panel.border = element_rect(linetype = "solid", fill = NA, colour = "black", size = 2),
                                        panel.background = element_blank(),               # element_rect(fill = NULL, colour = NULL),   # size = 0.5, linetype = "solid"
                                        panel.grid.major = element_blank(),               # element_line(colour = NULL),                # size = 0.5, linetype = 'solid',  
                                        panel.grid.minor = element_blank(),               # element_line(colour = NULL),                # size = 0.25, linetype = 'solid', 
                                        axis.ticks = element_line(colour = "black", size=2), 
                                        axis.ticks.length = unit(0.15, "cm"), 
                                        legend.position = "right",
                                        legend.key.size = unit(0.7, 'cm'),
                                        legend.title = element_text(size=20),
                                        legend.text = element_text(size=15)
                                    ) +
                                    geom_segment(aes(x = 1, y = 0, xend = seqTargetLength, yend = 0), size=1.0) +
                                    guides(fill=guide_legend(title="Strand")) +
                                    labs(x = "EVE sequence positions (bp)", y = "Coverage (Number of reads)")


# outfile
filename_Coverage_allREADS_PLOT <- paste0(output_path,"/",sample_name,".coverage.sRNA_sequence_all_reads.pdf");

pdf(file=filename_Coverage_allREADS_PLOT, width=24, height=10)
sample_cov_plot_both_strands
dev.off()





# WHOLE EVE INTEGRATION SEQUENCE: ONLY 21 SEQ-NTS SIZE
# -----------------------------------------------------------------------------

# extract mapped reads from specific sizes:
sRNA_sequence_mapped_reads_21nts.plot <- coverage_data2plot_updated %>% dplyr::filter(alignmentLength==21)
sRNA_sequence_mapped_reads_26_30nts.plot <- coverage_data2plot_updated %>% dplyr::filter(alignmentLength>=26 & alignmentLength<=30)


# coverage plot: 21nts 
sample_cov_plot_21nts <- ggplot(sRNA_sequence_mapped_reads_21nts.plot, aes(fill=Strand, y=Frequency_byStrand, x=Postition)) + 
                            geom_bar(position="dodge", stat="identity", , width=10.0) +
                            scale_fill_manual(values= c("red", "blue"),
                                            breaks= c("+", "-"),
                                            labels= c("+","-")
                            ) +
                            scale_x_continuous(expand = c(0, 0.1), limits = c(1,seqTargetLength), breaks = seq(0,seqTargetLength, by = 100)) +
                            scale_y_continuous(expand = c(0.005, 0), limits = c(-30,80), breaks = seq(-30,80, by = 10)) +
                            theme_minimal() +
                            theme(axis.text.x = element_text(angle = 0, hjust = 0.50, vjust = 1.0, colour = "black", size=20),
                                axis.text.y = element_text(colour = "black",size=20),
                                axis.title.x = element_text(color="black", size=23, face="bold"),
                                axis.title.y = element_text(color="black", size=23, face="bold"),
                                axis.text=element_text(size=10),
                                panel.border = element_rect(linetype = "solid", fill = NA, colour = "black", size = 2),
                                panel.background = element_blank(),               # element_rect(fill = NULL, colour = NULL),   # size = 0.5, linetype = "solid"
                                panel.grid.major = element_blank(),               # element_line(colour = NULL),                # size = 0.5, linetype = 'solid',  
                                panel.grid.minor = element_blank(),               # element_line(colour = NULL),                # size = 0.25, linetype = 'solid', 
                                axis.ticks = element_line(colour = "black", size=2), 
                                axis.ticks.length = unit(0.15, "cm"), 
                                legend.position = "right",
                                legend.key.size = unit(0.7, 'cm'),
                                legend.title = element_text(size=20),
                                legend.text = element_text(size=15)
                            ) +
                            geom_segment(aes(x = 1, y = 0, xend = seqTargetLength, yend = 0), size=1.0) +
                            guides(fill=guide_legend(title="Strand")) +
                            labs(x = "EVE sequence positions (bp)", y = "Coverage (Number of reads)")


# outfile
filename_Coverage_READS_21nts_PLOT <- paste0(output_path,"/",sample_name,".coverage.sRNA_sequence_21nts_reads.pdf");

pdf(file=filename_Coverage_READS_21nts_PLOT, width=24, height=10)
sample_cov_plot_21nts
dev.off()







# coverage plot: 26-30 nts 
sample_cov_plot_26_30nts <-  ggplot(sRNA_sequence_mapped_reads_26_30nts.plot, aes(fill=Strand, y=Frequency_byStrand, x=Postition)) + 
                                geom_bar(position="dodge", stat="identity", , width=10.0) +
                                scale_fill_manual(values= c("red", "blue"),
                                                breaks= c("+", "-"),
                                                labels= c("+","-")
                                ) +
                                scale_x_continuous(expand = c(0, 0.1), limits = c(1,seqTargetLength), breaks = seq(0,seqTargetLength, by = 100)) +
                                scale_y_continuous(expand = c(0.005, 0), limits = c(-30,80), breaks = seq(-30,80, by = 10)) +
                                theme_minimal() +
                                theme(axis.text.x = element_text(angle = 0, hjust = 0.50, vjust = 1.0, colour = "black", size=20),
                                    axis.text.y = element_text(colour = "black",size=20),
                                    axis.title.x = element_text(color="black", size=23, face="bold"),
                                    axis.title.y = element_text(color="black", size=23, face="bold"),
                                    axis.text=element_text(size=10),
                                    panel.border = element_rect(linetype = "solid", fill = NA, colour = "black", size = 2),
                                    panel.background = element_blank(),               # element_rect(fill = NULL, colour = NULL),   # size = 0.5, linetype = "solid"
                                    panel.grid.major = element_blank(),               # element_line(colour = NULL),                # size = 0.5, linetype = 'solid',  
                                    panel.grid.minor = element_blank(),               # element_line(colour = NULL),                # size = 0.25, linetype = 'solid', 
                                    axis.ticks = element_line(colour = "black", size=2), 
                                    axis.ticks.length = unit(0.15, "cm"), 
                                    legend.position = "right",
                                    legend.key.size = unit(0.7, 'cm'),
                                    legend.title = element_text(size=20),
                                    legend.text = element_text(size=15)
                                ) +
                                geom_segment(aes(x = 1, y = 0, xend = seqTargetLength, yend = 0), size=1.0) +
                                guides(fill=guide_legend(title="Strand")) +
                                labs(x = "EVE sequence positions (bp)", y = "Coverage (Number of reads)")


# outfile
filename_Coverage_READS_26_30nts_PLOT <- paste0(output_path,"/",sample_name,".coverage.sRNA_sequence_26-30nts_reads.pdf");

pdf(file=filename_Coverage_READS_26_30nts_PLOT, width=24, height=10)
sample_cov_plot_26_30nts
dev.off()








# WHOLE EVE SEQUENCE COVERAGE: ALL COVERAGE FIGURES MERGED
# -----------------------------------------------------------------------------
library(ggpubr)


title_plot <- paste0("EVE sequence coverage: ",sample_name);

sample_cov_plot_both_strands  <- sample_cov_plot_both_strands + xlab(NULL) + ylab("All mapped reads size");
sample_cov_plot_21nts         <- sample_cov_plot_21nts        + xlab(NULL) + ylab("Read size: 21 nts");
sample_cov_plot_26_30nts      <- sample_cov_plot_26_30nts     + xlab(NULL) + ylab("Read size: 26-30 nts");


figure <- ggarrange(sample_cov_plot_both_strands, 
                    sample_cov_plot_21nts,
                    sample_cov_plot_26_30nts,
                    ncol = 1, nrow = 3) 




# outfile
filename_all_coverages_combined <- paste0(output_path,"/",sample_name,".coverage.sRNA_sequence_all_coverage_plots.pdf");

pdf(file=filename_all_coverages_combined, width=24, height=13)

#figure
annotate_figure(figure,
                top = text_grob(title_plot, color = "black", face = "bold", size = 30),
                left = text_grob("Coverage (Reads Number)", color = "black", face = "bold", size = 25, rot = 90)
                )
dev.off()




# save this session
filename_Rsession <- paste0(output_path,"/Perdomo_etal_2025.small_RNAs.seq_size_distribution_and_coverage.",sample_name,".Rsession.RData");

save.image(file = filename_Rsession);




















#-@-#   ########################################################################################################################################################################
#-@-#   ## INFORMATION SAM/BAM FILES: SHORT INTRO TO FORMAT AND THEIR DATA STRUCTURE
#-@-#   ########################################################################################################################################################################
#-@-#
#-@-#
#-@-#   SAM/BAM file format:
#-@-#   
#-@-#   --------------------------------------------------------------------------------------------------
#-@-#   Col	Field	Type		Regexp/Range			Brief description
#-@-#   --------------------------------------------------------------------------------------------------
#-@-#   1	QNAME	String		[!-?A-~]{1,254}         	Query template NAME
#-@-#   2	FLAG	Int	    	[0, 216 − 1] 	        	bitwise FLAG
#-@-#   3	RNAME	String		\*|[:rname:∧*=][:rname:]* 	Reference sequence NAME12
#-@-#   4	POS	Int	    	[0, 231 − 1] 			1-based leftmost mapping POSition
#-@-#   5	MAPQ	Int	    	[0, 28 − 1] 			MAPping Quality
#-@-#   6	CIGAR	String		\*|([0-9]+[MIDNSHP=X])+ 	CIGAR string
#-@-#   7	RNEXT	String		\*|=|[:rname:∧*=][:rname:]* 	Reference name of the mate/next read
#-@-#   8	PNEXT	Int	    	[0, 231 − 1] 			Position of the mate/next read
#-@-#   9	TLEN	Int	    	[−231 + 1, 231 − 1] 		observed Template LENgth
#-@-#   10	SEQ	String		\*|[A-Za-z=.]+ 		segment SEQuence
#-@-#   11	QUAL	String		[!-~]+ 			ASCII of Phred-scaled base QUALity+33
#-@-#   --------------------------------------------------------------------------------------------------
#-@-#   
#-@-#   
#-@-#   TIP: Single end reads yield only 3 possible flag values: 0,4 and 16.
#-@-#    
#-@-#   * 0 means the read aligned in the forward direction: POSITIVE 
#-@-#   * 16 mean it aligned in the reverse direction: NEGATIVE 
#-@-#   * 4 means it didn't align.
#-@-#   
#-@-#   ########################################################################################################################################################################


