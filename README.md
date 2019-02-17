# Project BBT015

This project is a part of the course Advanced Bioinformatics (BBT015) at Chalmers. The aim with the project is to 
replicate a published bioinformatic paper. 

The project focues on the RNA-sequence part of a [paper](https://www.frontiersin.org/articles/10.3389/fmicb.2018.02929/)  with the name; *Transcriptional Landscape of a bla_KPC-2 Plasmid and Response 
to Imipenem Exposure in Escherichia coli TOP10* [1]. In the RNA-seq part of this study six *E.coli*-samples were investigated, 
and all these had been transformed with the pBIC-1a plasmid (conveys resicentence to antibiotics). Three samples acted as controlls 
while the other three had been exposed to the broad spectra Imipenem. The study investigated differences in expression levels 
between these two sample-groups, and overall it identified over 1500 differentially expressed RNA. Thus the goal with this 
project is try to replicate this, with other words the goal is to identify the same differentially RNA:s as the authors of the
paper did. 

# Project structure 

Each directory contains a readme.md file describing the role of that directory. The role of each directory
can be summarised as:

* *data*: Contains the sample and reference data. 
* *Docs*: Contains documentation. 
* *Scratch*: Contains files that can be saftley removed.
* *Results*: Contains the tables and figures produced by the code. 
* *Code*: Contains the Python, R and Shell-code required to replicate the project. 
* *Intermediate*: Intermediate files that aren't directly result. 

The source-code for the LaTex-writtenreport describing the result of the project can be found 
in the Docs folder. 


# References 

1. Jousset Agnès B, Rosinski-Chupin Isabelle, Takissian Julie, Glaser Philippe, Bonnin Rémy A, 
Naas Thierry. Transcriptional Landscape of a blaKPC-2 Plasmid and Response to Imipenem Exposure 
in Escherichia coli TOP10. *Frontiers in Microbiology*, 9(1):2929, dec 2018 
