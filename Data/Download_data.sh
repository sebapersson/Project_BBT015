#!/bin/bash

# This script will download the data required for the project and decompress it.
# More specifically this file will download:
# The 6 sample FASTQ-files.
# The 2 reference files (E.coli K12 and plasmid).
# The script won't download a file it is already present.
# It's important that the script is run from the data-directory!

# ------------------------- # Sample data # ---------------------------------- #
# Create Sample_data directory of not present
if [ ! -d "Sample_data" ]; then
    echo "Created Sample_data directory"
    mkdir Sample_data
fi

# Downloading the data into Sample_data
cd Sample_data/

# Function that will download a data-file if not present. If the file is compressed
# the function will decompress it.
# Input (order important)
# Arg1: Compressed file name
# Arg2: Download link
# Arg3: Name downloaded file
# Arg4: Sample number

download_sample_data ()
{
    # First check if sample is present 
    if [ -f $1 ]; then
	echo "Sample $4 already present"
    else
	echo "Sample $4 not present, will download"
	wget $2
	# Rename
	mv $3 ./"$1"
    fi

}

# Downloading the data
# Sample 1
dataLink1=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR278/001/ERR2780171/ERR2780171.fastq.gz
download_sample_data "Sample1.fastq.gz" $dataLink1 "ERR2780171.fastq.gz" 1

# Sample 2
dataLink2=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR278/002/ERR2780172/ERR2780172.fastq.gz
#download_sample_data "Sample2.fastq.gz" $dataLink2 "ERR2780172.fastq.gz" 2

# Sample 3
dataLink3=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR278/003/ERR2780173/ERR2780173.fastq.gz
#download_sample_data "Sample3.fastq.gz" $dataLink3 "ERR2780173.fastq.gz" 3

# Sample 4
dataLink4=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR278/004/ERR2780174/ERR2780174.fastq.gz
#download_sample_data "Sample4.fastq.gz" $dataLink2 "ERR2780174.fastq.gz" 4

# Sample 5
dataLink5=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR278/005/ERR2780175/ERR2780175.fastq.gz
#download_sample_data "Sample5.fastq.gz" $dataLink5 "ERR2780175.fastq.gz" 5

# Sample 6
dataLink6=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR278/006/ERR2780176/ERR2780176.fastq.gz
#download_sample_data "Sample2.fastq.gz" $dataLink2 "ERR2780176.fastq.gz" 6

# Move back to data directory 
echo "All samples are downloaded"
cd ..


