#!/bin/bash

# Move to Intermediate directory
cd ../../Intermediate/
echo ""
# Create Alignment data directory and subdirectories if non-existent 
if [ ! -d "Alignment_data" ]; then
    echo "Creating Alignment_data directory"
    mkdir Alignment_data
    cd Alignment_data
    echo "Creating Alignment_data subdirectories E_coli and Plasmid" 
    mkdir E_coli
    mkdir Plasmid
    echo "Directories created!"
else
    cd Alignment_data
    echo "Alignment data directory already exists"
fi
echo ""

# Function for Aligning data 
# Variables:
# $1 is sample file
# $2 unzipped sample file
# $3 is samfile for E_coli genome alignment
# $4 is bamfile for E_coli genome alignment
# $5 is samfile for plasmid alignment
# $6 is bamfile for plasmid alignment
# $7 is sample number
data_alignment_function()
{
    # Perform alignments for E.coli genome and plasmid if plasmid bamfile does not exist
    cd E_coli
    if [ ! -f  $4 ]; then
	cd ..

	echo "Retrieving Sample $7" 
	#Retrieve data
	cp ../../Data/Sample_data/$1 .
	echo "Sample $7 retrieved!"
	echo ""
	echo "Unzipping Sample $7"
	gunzip $1
	chmod 777 $2
	echo "Sample $7 unzipped!"
	echo ""
	
	#Create alignment for E.coli genome 
	cd E_coli
	echo "Creating alignment of Sample $7 to E.coli genome. Takes about 10 minutes"
	./../../../Bin/bowtie-0.12.7/bowtie -tS E_coli_index ../$2 $3
	echo "Alignment created!"
	echo ""
	echo "Converting to bam file format" 
	samtools view -Sb $3 > $4
	echo "Alignment file is now a bam file!"

	#Remove sam file to clear up space 
	rm $3
	
	#Create alignment for plasmid
	cd ../Plasmid
	echo "Creating alignment of Sample $7 to plasmid. Takes about 10 minutes"
	./../../../Bin/bowtie-0.12.7/bowtie -tS pBIC_1a_index ../$2 $5
	echo "Alignment created!"
	echo ""
	echo "Converting to bam file format" 
	samtools view -Sb $5 > $6
	echo "Alignment file is now a bam file!"

	#Remove sam file to clear up space 
	rm $5
	
	#Remove unzipped sample file
	cd ..
	echo "Removing unzipped file of Sample $7"
	rm $2
	echo "File removed!"
	echo ""
	echo "Alignment performed successfully for Sample $7"
    else
	echo "Alignment files for Sample $7 already exists"
	# Move back to Alignment_data directory
	cd ..
    fi
}

# Sample 1
data_alignment_function "Sample1.fastq.gz" "Sample1.fastq" "Sample1.map.sam" "Sample1.map.bam" "Sample1_plasmid.map.sam" "Sample1_plasmid.map.bam" 1 

# Sample 2
data_alignment_function "Sample2.fastq.gz" "Sample2.fastq" "Sample2.map.sam" "Sample2.map.bam" "Sample2_plasmid.map.sam" "Sample2_plasmid.map.bam" 2 

# Sample 3
data_alignment_function "Sample3.fastq.gz" "Sample3.fastq" "Sample3.map.sam" "Sample3.map.bam" "Sample3_plasmid.map.sam" "Sample3_plasmid.map.bam" 3 

# Sample 4
data_alignment_function "Sample4.fastq.gz" "Sample4.fastq" "Sample4.map.sam" "Sample4.map.bam" "Sample4_plasmid.map.sam" "Sample4_plasmid.map.bam" 4 

# Sample 5
data_alignment_function "Sample5.fastq.gz" "Sample5.fastq" "Sample5.map.sam" "Sample5.map.bam" "Sample5_plasmid.map.sam" "Sample5_plasmid.map.bam" 5 

# Sample 6
data_alignment_function "Sample6.fastq.gz" "Sample6.fastq" "Sample6.map.sam" "Sample6.map.bam" "Sample6_plasmid.map.sam" "Sample6_plasmid.map.bam" 6 

exit 0



