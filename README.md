# Chrom-Lasso
Chrom-Lasso is a tool to analyze Hi-C data for identifying chromatin interactions.
The code of Chrom-Lasso is made up of Cpp and R.
The following files are code, input files, test data, and tutorial for Chrom-Lasso.
1. Code
This folder contains code for Chrom-Lasso, the usage of them is shown in the "Tutorial" folder.
2. Prepare_Input_File
This folder contains input files needed by Chrom-Lasso to do analysis, mainly are domain files and cutting site files for mouse and human,
The domain file are constant for specific species, but the genomic coordinate in the file can be changed with genome version by tool LiftOver.
The cutting site file can be produced by tool Oligomatch with reference genome and cutting sequence for the restriction endonuclease.
3. Test_Data
This folder contains "sortChr" files for mouse and human for testing Chrom-Lasso, the "sortChr" file are preprocessed file from raw "fastq" sequencing files,
the preprocessing of "fastq" files can be done with our totorial, which is the same as JUICER, so you can use JUICER as well.
The JUICER generates "merged_nodups" file, and you can use the Shell scripts in totorial that sorts the file by chromosome order to generate "sortChr" file for further analysis.
4. Tutorial
This folder contains the analysis pipeline for Mouse and Human, you can directly use "sotrChr" files in "Test_Data" folder to test the pipeline, 
you can also use your own "fastq" files from Hi-C sequencing to test the totorial step by step from the preprocessing to running Chrom-Lasso.

 
