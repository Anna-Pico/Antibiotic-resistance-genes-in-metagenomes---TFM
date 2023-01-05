# Antibiotic resistance genes in metagenomes - TFM

This repository contains the shell scripts and the database I used to analyze wastewater metagenomes and identify their resistome. The resistome is the set of antibiotic resistance genes (ARG) present within a sample. the basic scheme of the proces I followed is the following:

![image](https://user-images.githubusercontent.com/116899896/203780212-5de03fc8-f51a-4d13-b3f1-62d22ceece3b.png)

Tools used to achieve this goals are 

- Fastx-toolkit (Hannon, G.J., 2010): https://github.com/agordon/fastx_toolkit
- Diamond (Buchfink, B., et al., 2021): https://github.com/bbuchfink/diamond
- Metaxa (Bengtsson, J., et al., 2011): https://microbiology.se/software/metaxa2/

And it is structured in 3 files:

1) *analisi de sequencies.sh*: this shell script contains the code for the analysis of one metagenome. 
2) *referenceARG.dmnd*: The database that contains the reference sequences for the antibiotic resistance genes. It is based in the ARGminer-v1.1.1-A (2019) and it has been slightly modified to facilitate the agrupation by class of antibiotic (Arango-Argoty, G. A., et al., 2020).
3) *loop.sh* This file allows to process continuously all the fastq files present in the folder.

The required input is a fastq file containing all the genomic sequences (metagenome). The output is a csv file with two columns that contains the counts of reads related to each antibiotic class. Besides, we obtain raw outputs from the Diamond and Metaxa tools. 

\* This repositore also includes an R code that has been used to analyse 56 resistomes from wastewater samples. Although the code is data-dependent and can not be run by a third parties it may be useful to better understand our work or can be an inspiration for future studies.

References:

- Arango-Argoty, G. A., Guron, G. K. P., Garner, E., Riquelme, M. V, Heath, L. S., Pruden, A., Vikesland, P. J., & Zhang, L. (2020). ARGminer: a web platform for the crowdsourcing-based curation of antibiotic resistance genes. Bioinformatics, 36(9), 2966–2973. https://doi.org/10.1093/bioinformatics/btaa095
- Bengtsson, J., Eriksson, K. M., Hartmann, M., Wang, Z., Shenoy, B. D., Grelet, G. A., Abarenkov, K., Petri, A., Alm Rosenblad, M., & Nilsson, R. H. (2011). Metaxa: A software tool for automated detection and discrimination among ribosomal small subunit (12S/16S/18S) sequences of archaea, bacteria, eukaryotes, mitochondria, and chloroplasts in metagenomes and environmental sequencing datasets. Antonie van Leeuwenhoek, International Journal of General and Molecular Microbiology, 100(3), 471–475. https://doi.org/10.1007/s10482-011-9598-6
- Buchfink, B., Reuter, K., & Drost, H.-G. (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. Nature Methods. https://doi.org/10.1038/s41592-021-01101-x
- Hannon, G. J. (2010). FASTX-Toolkit. http://hannonlab.cshl.edu/fastx_toolkit


