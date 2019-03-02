# About

This directory contains all the code required to produce the result from the RNA-seq data in the *data* directory. In the following sections is a brief explanation of each file present in this directory. 

## Shell-code

### Run_all.sh

This script will run the entire analysis and produce all tables and figures, with other words it will fully replicate the data-analysis. Thus to fully replicate the project one only has to clone this repository from GitHub and run the run-all script. The results from will script is stored in Results/Figures_copy and Results/Tables_copy.

The different steps performed by this scripts are:

1. Downloading the sample data (and storing in the *Data*-directory).
2. Creating an extra annotation file for the plasmid. 
3. Downloading bowtie-0.12.7 to the *Bin*-directory.
4. Aligning the reads to and storing the result in Intermediate/Count_mat using bowtie. 
5. Checking that the required R-packages are installed. 
6. Creating count matrices for the plasmid and *E. coli* and storing them in *Intermediate*-directory.
7. Performing DESeq analysis. 

### Download_data.sh

This script will download the sample data required for the project and store it in the *Data*-directory. More specifically this file will download the 6 sample FASTQ-files (if they already aren't present).

### Download_bowtie-0_12_7.sh

This script will download bowtie version 0.12.7 to the *Bin* directory if it's not present already. If bowtie is downloaded due to not being present the script will check the installation. Besides downloading bowtie the script will create the bowtie indices for the *E. coli* and the plasmid references. 

### Alignment.sh

This script will align the sample data (present in the *Data/Sample_data*-directory) to the reference genomes. The alignment takes approximately 10 minutes per sample. Altogether 12 bam-files will be created (six for the *E. coli* and six for the plasmid) and stored in the Intermediate/Alignment_data*. If the bam-files already are present the script won't re-align the data.

## R-code 

### Differential_analysis_E_coli.R and Differential_analysis_plasmid.R

These files performs the differential analysis of for the plasmid and *E .coli* using the *DEseq2*-package. Output from these files are:

* Table with most significant genes (adjusted false discovery rate p-value < 0.05).
* Heat-map of difference between samples (based on Poisson distance).
* PCA-plot for the different samples, the projected data normalized by the vst-method.
* Volcano plot over all genes.
* A histogram over the p-values from the Wald-test performed by DESeq2. 

These figures and tables can be found in *Results/Tables_copy* and *Results/Figures_copy*. 

### Create_count_mat_plasmid.R and Create_count_mat_E_coli.R

These file will create the count matrix for *E.coli* and the plasmid from the aligned bam-files located in the *Intermediate/Alignment_data/E_coli* folder. The resulting count matrices are stored in the *Intermediate/Count_mat*-directories (which will be created if they aren't present. 

### Check_libraries_installed.R

This program checks that all the R-packages required for the differential analysis are installed. If a package isn't installed the program will install it locally. Since dependencies are tricky among the *bioconductor* packages (e.g *DESeq2*) the program will exit if a *bioconductor* isn't installed (and it won't try to installing it).

## Python-code 

### Annotate_plasmid.py

This file extracts the Gene-name and product (gene description) from the *pBIC_1a.gff3* file in the data-directory using regular expression. The result is then written the results to a tab-separated file in the *Data/Reference_data/Plasmid folder*. The annotation file is used to annotate (add gene description) for the plasmid-genes in the differential analysis step. 



