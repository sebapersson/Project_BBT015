#!/bin/bash

# Move to Bin directory
cd ../../Bin

# If bowtie directory does not exist, download bowtie-0.12.7 and unzip
# When bowtie zip file is unzipped, a directory called bowtie-0.12.7 is created
# Zip file is removed
if [ ! -d "bowtie-0.12.7" ]; then
    echo "Downloading bowtie-0.12.7"    
    wget --content-disposition -c https://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip/download &> /dev/null
    unzip bowtie-0.12.7-linux-x86_64.zip > /dev/null
    rm bowtie-0.12.7-linux-x86_64.zip
    echo "Download complete!"
else
    echo "bowtie-0.12.7 already present"
fi


# Move to bowtie indexes folder to create reference genome indexes
cd bowtie-0.12.7/indexes

# If E.coli index does not exist, create it
if [ ! -f "E_coli_index.1.ebwt" ]; then
    echo "Creating index for E_coli genome"
    ./../bowtie-build ../../../Data/Reference_data/E_coli/E_coli.fasta E_coli_index > /dev/null
    echo "Index built!"
else
    echo "E_coli genome index already present"
fi

# If plasmid index does not exist, create index
if [ ! -f "pBIC_1a_index.1.ebwt" ]; then
    echo "Creating index for plasmid pBIC_1a" 
    ./../bowtie-build ../../../Data/Reference_data/Plasmid/pBIC_1a.fasta pBIC_1a_index > /dev/null
    echo "Index built!"
else
    echo "Index for plasmid pBIC_1a already present" 
fi

echo ""
echo ""

# Test installment by alignment of random sequence
echo "Testing alignment of random sequence to E.coli index"
./../bowtie -c E_coli_index GCGTGAGCTATGAGAAAGCGCCACGCTTCC 

echo ""

# Test installment by alignment of random sequence
echo "Testing alignment of random sequence to plasmid index"
./../bowtie -c pBIC_1a_index GCGTGAGCTATGAGAAAGCGCCACGCTTCC 

echo ""

if [ -f "pBIC_1a_index.1.ebwt" ]; then
    echo "Download of bowtie-0.12.7 and indexes complete!"
else
    echo "Something went wrong"
fi
