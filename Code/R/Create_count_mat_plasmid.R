# This file will create the count matrix for the plamis from the aligned bam-files located in the 
# Intermediate/Alignment_data/Plasmid folder 
# Note the script should be run from the Code/R directory for all file-paths to correct. 
rm(list = ls())

# Quit if countMatrix already exists 
condVec <- logical(3)
condVec[1] <- dir.exists("../../Intermediate/Count_mat/Plasmid")
condVec[2] <- file.exists("../../Intermediate/Count_mat/Plasmid/assays.h5")
condVec[3] <- file.exists("../../Intermediate/Count_mat/Plasmid/se.rds")

if(all(condVec)){
  print("Count matrix for E. coli is already present, exit status 0")
  quit(status = 0)
}

library(GenomicAlignments)
library(Rsamtools)
library(GenomicFeatures)
library(HDF5Array)

# Creating the file-paths
bamFilesDir <- "../../Intermediate/Alignment_data/Plasmid/"
bamFilePaths <- c(paste0(bamFilesDir, "Sample1_plasmid.map.bam"), paste0(bamFilesDir, "Sample2_plasmid.map.bam"), 
                  paste0(bamFilesDir, "Sample3_plasmid.map.bam"), paste0(bamFilesDir, "Sample4_plasmid.map.bam"), 
                  paste0(bamFilesDir, "Sample5_plasmid.map.bam"), paste0(bamFilesDir, "Sample6_plasmid.map.bam"))

# Check that the file exists                
if(!all(file.exists(bamFilePaths))){
  print("Error with bamFilePaths, program will exit.")
  quit(status = 1)
}

# Creating a reference to the BAM-files
# Depending on computer yield-size (how many reads processed at a time) might need to be changed
bamFiles <- BamFileList(bamFilePaths, yieldSize=2000000)

# Annotation data 
filePathGFF <- "../../Data/Reference_data/Plasmid/pBIC_1a.gff3"

if(!file.exists(filePathGFF)){
  print("Error with GFF file-path E.coli, program will exit.")
  quit(status = 1)
}

# Creating TxDb object of annoation data, by default assumes circular genome
annotationData <- makeTxDbFromGFF(filePathGFF, format = "gff3")

# Creating GRangesList of all exons grouped by genes 
geneList <- exonsBy(annotationData, by="gene")

# Creating the count Matrix for E.coli
# Will create counts by genes
# Use the most convservitaive union
# Data is not paired end 
countMatPlasmid <- summarizeOverlaps(features=geneList, reads=bamFiles,
                                   mode="Union",
                                   singleEnd=TRUE,
                                   ignore.strand=TRUE,
                                   fragments=FALSE )

# Saving the count-matrix, create directoires if they don't exist 
if(!dir.exists("../../Intermediate/Count_mat")){
  dir.create("../../Intermediate/Count_mat")
}

dirToSaveCountMat <- "../../Intermediate/Count_mat/Plasmid/"
saveHDF5SummarizedExperiment(countMatPlasmid, dir=dirToSaveCountMat, replace=TRUE)

print("The count matrix for E.coli is now located in Intermediate/Count_mat/E_coli")

quit(status = 0)

