# tRNA_reference_construction
You will find here the different pipelines used by the NGS Core Facility "EpiRNA-Seq" from Nancy, France (https://umsibslor.univ-lorraine.fr/en/facility/epitranscriptomics-sequencing-epirna-seq) to create tRNA references from genomic tRNAs predictions.

!!! Be aware !!! All the scripts has been used in a local computer, with relatives paths for the inputs and outputs. If you want to use these scripts on your own, please be careful of the inputs and directories names along the script.

Every manual input that you should modify will be indicated by this symbol: //!\\\

# Requirements

An Unix shell and the R software are needed to use these scripts.

The "Unix_scripts" needs Trimmomatic, Bowtie2 and samtools. You can naturally use your own trimmming and alignment tools if you prefer, but be sure that the trimming output is in FASTQ format, and the alignment output in SAM format.

# Pipeline workflow
