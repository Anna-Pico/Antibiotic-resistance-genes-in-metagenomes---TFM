# Antibiotic resistance genes in metagenomes - TFM

This repository contains the shell scripts and the database I used to analyze wastewater metagenomes and identify their resistome. The resistome is the set of antibiotic resistance genes (ARG) present within a sample. the basic scheme of the proces I followed is the following:

![image](https://user-images.githubusercontent.com/116899896/203780212-5de03fc8-f51a-4d13-b3f1-62d22ceece3b.png)

Tools used to achieve this goals are 

- Fastx-toolkit https://github.com/agordon/fastx_toolkit
- Diamond https://github.com/bbuchfink/diamond
- Metaxa https://microbiology.se/software/metaxa2/

And it is structured in 3 files:

1) *analisi de sequencies.sh*: this shell script contains the code for the analysis of one metagenome. 
2) *referenceARG.dmnd*: The database that contains the reference sequences for the antibiotic resistance genes. It is based in the ARGminer-v1.1.1-A (2019) and it has been slightly modified to facilitate the agrupation by class of antibiotic.
3) *loop.sh* This file allows to process continuously all the fastq files present in the folder.

The required input is a fastq file containing all the genomic sequences (metagenome). The output is a csv file with two columns that contains the counts of reads related to each antibiotic class. Besides, we obtain raw outputs from the Diamond and Metaxa tools. 
