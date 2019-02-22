#!/bin/bash

# This script will download the sample data required for the project and decompress it.
# More specifically this file will download:
# The 6 sample FASTQ-files.
# The script won't download a file it is already present.
# It's important that the script should run from the data-directory!
# The script will fix file-permission for all files

# ------------------------- # Sample data # ---------------------------------- #
# Create Sample_data directory of not present
if [ ! -d "Sample_data" ]; then
    echo "Created Sample_data directory"
    mkdir Sample_data
fi

# Downloading the data into Sample_data
cd Sample_data/

# Function that will download a data-file if not present. If not it will 
# be downloaded. 
# Input (order important)
# Arg1: Compressed file name
# Arg2: Download link
# Arg3: Name downloaded file
# Arg4: Sample number

download_sample_data ()
{
    if [ -f $1 ]; then
	echo "Sample $4 already present"
    else
	echo "Sample $4 not present. Starting to download it"
	wget -q $2
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
download_sample_data "Sample2.fastq.gz" $dataLink2 "ERR2780172.fastq.gz" 2

# Sample 3
dataLink3=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR278/003/ERR2780173/ERR2780173.fastq.gz
download_sample_data "Sample3.fastq.gz" $dataLink3 "ERR2780173.fastq.gz" 3

# Sample 4
dataLink4=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR278/004/ERR2780174/ERR2780174.fastq.gz
download_sample_data "Sample4.fastq.gz" $dataLink4 "ERR2780174.fastq.gz" 4

# Sample 5
dataLink5=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR278/005/ERR2780175/ERR2780175.fastq.gz
download_sample_data "Sample5.fastq.gz" $dataLink5 "ERR2780175.fastq.gz" 5

# Sample 6
dataLink6=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR278/006/ERR2780176/ERR2780176.fastq.gz
download_sample_data "Sample6.fastq.gz" $dataLink6 "ERR2780176.fastq.gz" 6

# Data should be read only and non-compressed 
chmod 444 Sample1.fastq.gz && gunzip Sample1.fastq.gz
chmod 444 Sample2.fastq.gz && gunzip Sample2.fastq.gz
chmod 444 Sample3.fastq.gz && gunzip Sample3.fastq.gz
chmod 444 Sample4.fastq.gz && gunzip Sample4.fastq.gz
chmod 444 Sample5.fastq.gz && gunzip Sample5.fastq.gz
chmod 444 Sample6.fastq.gz && gunzip Sample6.fastq.gz

# Move back to data directory 
echo "All samples are downloaded"
cd ..


# --------------------------- # Reference data # ------------------------------ #
# Check if Reference_data is available, if not something is wrong!
if [ ! -d "Reference_data" ]; then
    echo "Reference data isn't present, check directory structure or download it."
    exit 1
fi

# If reference present, fix file-permissions, data should be read only
cd Reference_data/E_coli

# E.coli directory 
for file in $( ls ); do
    chmod 444 $file
done 

cd ../Plasmid

# Plasmid directory 
for file in $( ls ); do
    chmod 444 $file
done

exit 0
