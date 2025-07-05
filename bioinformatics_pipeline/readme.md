## Description of the pipeline.

The pipeline was applied to all mosquito samples described in the paper. Raw reads are publicly available in the NCBI SRA database under the ID SRAXXXXX.
Reads were mapped against two different reference sequences: 1) Cell Fusion Agent Virus genome RioPiedrasXXX, and 2) Endogenous Viral Element 3 from Crava etal., 2020.
See methods in the paper.

A short simple pipeline was created to combine all needed tasks: identify, measure sizes, quantify and locate the distribution of small RNAs of different lengths in the reference sequences separately.
To run properly the pipeline `main_pipeline.sh` script path files and programs/packages must be previously installed. Additional modifications in the R script
For each single sample, the command line 
