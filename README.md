# Project BBT015
## Description

This project is a part of the course Advanced Bioinformatics (BBT015) at Chalmers University of Technology. The aim with the project is to replicate a published bioinformatic paper. 

The project focues on the RNA-sequence part of a [paper](https://www.frontiersin.org/articles/10.3389/fmicb.2018.02929/)  with the name; *Transcriptional Landscape of a bla_KPC-2 Plasmid and Response to Imipenem Exposure in Escherichia coli TOP10* [1]. In the RNA-seq part of this study six *E.coli*-samples were investigated, 
and all these had been transformed with the pBIC-1a plasmid (conveys resisance to antibiotics). Three samples acted as controls 
whilst the other three had been exposed to the broad spectra Imipenem. The study investigated differences in expression levels 
between these two sample-groups, and overall it identified over 1500 differentially expressed RNA. Thus, the goal with this 
project is try to replicate this, with other words, to identify the same differentially RNA:s as the authors of the
paper did. 

To fully replicate the produced result in this project the requirements in the *Requirements for replication of result* section below should be fulfilled. Given this all the tables and figures in the *Results*-directory can be re-created by running *Run_all.sh* script in the *Code/Shell*-directory from its location. 

## Project structure 

Each directory contains a readme.md file describing the role of that directory. The role of each directory
can be summarised as:

* *Data*: Contains the sample and reference data. 
* *Docs*: Contains documentation. 
* *Scratch*: Contains files that can be saftley removed.
* *Results*: Contains the tables and figures produced by the code. 
* *Code*: Contains the Python, R and Shell-code required to replicate the project. 
* *Intermediate*: Intermediate files that aren't directly result. 
* *Bin*: Contains external scripts and compiled programs required to run the analysis. 

The source-code for the LaTeX-written report describing the result of the project can be found 
in the *Docs* folder. 

## Requirements for replication of result

This entire repository was created on Ubuntu Linux, and the code should be able to run on any Unix-based system. R version 3.5.2 and Python version 3.6.7 was used to produce the result. Since the python-code only relies on standard libraries (**sys**, **re** and **os**) a different Python-version shouldn't produce a different result than that in the *Results*-directory. A different R-version however might produce a different result, which is mainly due to the *bioconductor*-packages. Thus to ensure a full replication R version 3.5.2 and the package-versions below should be used. The non-standard R-packages used for the analysis were:

* **Rsamtools**, version 1.34.1.
* **GenomicAlignments**, version 1.18.1.
* **GenomicFeatures**, version 1.34.3
* **DESeq2**, version 1.22.2
* **pheatmap**, version 1.0.12
* **ggplot2**, version 3.1.0
* **PoiClaClu**, version 1.0.2.1

Besides these packages **xtable** (1.8.3) and **RColorBrewer** (1.1.2) were used. However different versions of those shouldn't affect the analysis in any way, since they are only used to produce tables for LaTex and colors for the graphs.

## References 

1. Jousset Agnès B, Rosinski-Chupin Isabelle, Takissian Julie, Glaser Philippe, Bonnin Rémy A, 
Naas Thierry. Transcriptional Landscape of a blaKPC-2 Plasmid and Response to Imipenem Exposure 
in Escherichia coli TOP10. *Frontiers in Microbiology*, 9(1):2929, dec 2018 
