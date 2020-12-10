# tRNA_reference_construction
You will find here the different pipelines used by the NGS Core Facility "EpiRNA-Seq" from Nancy, France (https://umsibslor.univ-lorraine.fr/en/facility/epitranscriptomics-sequencing-epirna-seq) to create tRNA references from genomic tRNAs predictions.

!!! Be aware !!! All the scripts has been used with relatives paths for the inputs and outputs. If you want to use these scripts on your own, please be careful of the inputs and directories names along the script.

Every manual input that you should modify will be indicated by this symbol: //!\\

# Requirements

An Unix shell and the R software are needed to use these scripts.

The "Unix_scripts" needs Trimmomatic, Bowtie2 and samtools. You can naturally use your own trimmming and alignment tools if you prefer, but be sure that the trimming output is in FASTQ format, and the alignment output in SAM format.

# How to create references (R script)
  1. Change the path of "ExtDataDir" to your current directory where you have a FASTA file with tRNA geneomic predictions.
  2. Check if "fasta_file_RNApattern" is correct
  3. (Optional) Change the "ResultsDir" path to store the output outside from the "extDataDir" by default
  4. Do NOT run everything! This script is written is 3 parts, and each should be run one after the other. 
  5. Done!

# How to validate references (Unix script) (for advanced users)
  1. Put the script inside the parent folders of your tRNA samples
  
         \__\ Sample_1 
         \__\ Sample_2
         \__\ Sample_3
         |_| Script HERE
  
  2. Change the "Genome" path to your different rRNA/tRNA references if needed ('/illumina/runs/iGenomes' by default)
  3. Change the paths to Bowtie2 and Trimmomatic
  4. The command line should be executed as follow :
      - For 1TA1 script : 'bash 1TA1_Validation_tRNA_reference <rRNA reference> <Non duplicated tRNA reference (Step 1)> <Merged tRNA reference (Step 2)>'
      - For 1TA2 script : 'bash 1TA2_Validation_tRNA_reference <Non duplicated tRNA reference (Step 1)> <Optimized tRNA reference (Step 3)>'
  5. Done!
  
  
If you need help or additional informations, please send a mail to fpichot@uni-mainz.de.
