# This file will create the count matrix for the plamis from the aligned bam-files located in the 
# Intermediate/Alignment_data/Plasmid folder 
# Note the script should be run from the Code/R directory for all file-paths to correct. 
rm(list = ls())

# Quit if countMatrix already exists 
condVec <- logical(3)
condVec[1] <- dir.exists("../../Intermediate/Count_mat/Plasmid")
condVec[2] <- file.exists("../../Intermediate/Count_mat/Plasmid/Count_mat_plasmid.dat")
condVec[3] <- file.exists("../../Intermediate/Count_mat/Plasmid/Sample_data_plasmid.dat")

if(all(condVec)){
  print("Count matrix and sample data for plasmid is already present")
  quit(status = 0)
}

library(GenomicAlignments)
library(Rsamtools)
library(GenomicFeatures)

# ----------------------------------- # Creating count matrix # --------------------------------------- #
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

# Creating the count Matrix for the plasmid 
# Will create counts by genes
# Use the most convservitaive union
# Data is not paired end 
countMatPlasmid <- summarizeOverlaps(features=geneList, reads=bamFiles,
                                   mode="Union",
                                   singleEnd=TRUE,
                                   ignore.strand=TRUE,
                                   fragments=FALSE )


# Reading the sample data
pathToSampleData <- "../../Data/Sample_data/Sample_data_info.dat"
if(!file.exists(pathToSampleData)){
  message("Sample data for E.coli doesn't exist.")
  quit(status = 1)
}

sampleData <- read.table(file = pathToSampleData, header = T, sep = "\t")

# Removing spaces for names
sampleData$Source.Name <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6")

# Filtering out names, condition and dose
subSampleData <- sampleData[, c(1, 28, 29)]
names(subSampleData) <- c("sample", "condition", "dose")
subSampleData$sample <- as.factor(subSampleData$sample)

# Matching row-names to sample names 
rownames(subSampleData) <- c("Sample1_plasmid.map.bam", "Sample2_plasmid.map.bam", "Sample3_plasmid.map.bam", 
                             "Sample4_plasmid.map.bam", "Sample5_plasmid.map.bam", "Sample6_plasmid.map.bam")
colData(countMatPlasmid) <- DataFrame(subSampleData)

# ---------------------------------- # Write count matrix to disk # ----------------------------------- #
# Saving the count-matrix, create directoires if they don't exist 
if(!dir.exists("../../Intermediate/Count_mat")){
  dir.create("../../Intermediate/Count_mat")
}

if(!dir.exists("../../Intermediate/Count_mat/Plasmid")){
  dir.create("../../Intermediate/Count_mat/Plasmid")
}

# Saving the count matrix and the colData
sampleInfo <- colData(countMatPlasmid)
countMatrix <- assay(countMatPlasmid)

write.table(sampleInfo, file = "../../Intermediate/Count_mat/Plasmid/Sample_data_plasmid.dat")
write.table(countMatrix, file = "../../Intermediate/Count_mat/Plasmid/Count_mat_plasmid.dat")
print("Count matrix created for plasmid")

quit(status = 0)

