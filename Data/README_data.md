# Data information 

For downloading the sample data, and store it in a new directory called *Sample_data*, run the script *Download_data.sh* located in the *Data* directory. A description of the different data-files is given below:

## Sample data 
The sample data consists of six RNA-seq samples from E.coli K-12 transformed with the pBIC-1a plasmid. Sequencing of the data was performed using Illumina HiSeq 2500 in singe-read mode with 50 cycles (not paried-ends). 

The first three samples (sample 1-3) act as control and have not been exposed to the broad spectra imipenem. The last three samples (sample 4-6) have been exposed to imipenem for 10 minutes before extraction of RNA. 

The sample data was taken from ArrayExpress, association number E-MTAB-7190, 2019-02-17. The links used for downloading the data are:

* Sample 1: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR278/001/ERR2780171/ERR2780171.fastq.gz
* Sample 2: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR278/002/ERR2780172/ERR2780172.fastq.gz
* Sample 3: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR278/003/ERR2780173/ERR2780173.fastq.gz
* Sample 4: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR278/004/ERR2780174/ERR2780174.fastq.gz
* Sample 5: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR278/005/ERR2780175/ERR2780175.fastq.gz
* Sample 6: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR278/006/ERR2780176/ERR2780176.fastq.gz

## Reference data 

The reference data consists of two parts, a reference genome for *E.coli* and a reference genome for the for the pBIC-1a plasmid. 

### *E.coli* reference
The strain used in the study is; *Escherichia coli* str. K-12 substr. DH10B. The reference genome for this strain was downloaded from [https://www.ncbi.nlm.nih.gov/nucleotide](NCBI) nucleotide data-base, accession number; CP000948.1. The downloaded files are:

* Sequence, fasta: Original file renamed to *E_coli.fasta*, downloaded 2019-02-20
* Genome annotation, GFF3: Original file renamed to *E_coli.gff3*, downloaded 2019-02-20. 

The files can be found in *Data/Reference_samples/E_coli* folder. 

### pBIC-a plasmid reference
The plasmid used in the study is; *Klebsiella pneumoniae* strain BIC-1 plasmid pBIC-1a. The reference genome for this plasmid was downloaded from [https://www.ncbi.nlm.nih.gov/nucleotide](NCBI) nucleotide data-base, accession number; CP022574. The downloaded files are:

* Sequence, fasta: Original file renamed to pBIC_1a.fasta, downloaded 2019-02-20
* Annotation, GFF3: Original file renamed to pBIC_1a.gff3, download 2019-02-20

The files can be found in *Data/Reference_samples/Plasmid* folder. 

