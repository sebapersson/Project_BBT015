# This file will create the count matrix form E.coli from the aligned bam-files located in the 
# Intermediate/Alignment_data/E_coli folder 
# The script will save the resulting Deseq object. 
# Note the script should be run from the Code/R directory for all file-paths to correct. 
rm(list = ls())

condVec <- logical(3)
condVec[1] <- dir.exists("../../Intermediate/Count_mat/E_coli")
condVec[2] <- file.exists("../../Intermediate/Count_mat/E_coli/Count_mat_E_coli.dat")
condVec[3] <- file.exists("../../Intermediate/Count_mat/E_coli/Sample_data_E_coli.dat")

if(all(condVec)){
  print("Count matrix for E. coli is already present, exit status 0")
  quit(status = 0)
}

library(GenomicAlignments)
library(Rsamtools)
library(GenomicFeatures)

# ----------------------------------- # Creating count matrix # --------------------------------------- #
# Creating the file-paths
bamFilesDir <- "../../Intermediate/Alignment_data/E_coli/"
bamFilePaths <- c(paste0(bamFilesDir, "Sample1.map.bam"), paste0(bamFilesDir, "Sample2.map.bam"), 
                  paste0(bamFilesDir, "Sample3.map.bam"), paste0(bamFilesDir, "Sample4.map.bam"), 
                  paste0(bamFilesDir, "Sample5.map.bam"), paste0(bamFilesDir, "Sample6.map.bam"))

# Check that the file exists                
if(!all(file.exists(bamFilePaths))){
  print("Error with bamFilePaths, program will exit.")
  quit(status = 1)
}

# Creating a reference to the BAM-files
# Depending on computer yield-size (how many reads processed at a time) might need to be changed
bamFiles <- BamFileList(bamFilePaths, yieldSize=2000000)
               
# Annotation data 
filePathGFF <- "../../Data/Reference_data/E_coli/E_coli.gff3"

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
countMatEcoli <- summarizeOverlaps(features=geneList, reads=bamFiles,
                                   mode="Union",
                                   singleEnd=TRUE,
                                   ignore.strand=TRUE,
                                   fragments=FALSE )

# ----------------------------------- # Annotate count matrix # --------------------------------------- #
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
rownames(subSampleData) <- c("Sample1.map.bam", "Sample2.map.bam", "Sample3.map.bam", 
                             "Sample4.map.bam", "Sample5.map.bam", "Sample6.map.bam")
colData(countMatEcoli) <- DataFrame(subSampleData)


# ---------------------------------- # Write count matrix to disk # ----------------------------------- #
# Saving the count-matrix, create directoires if they don't exist 
if(!dir.exists("../../Intermediate/Count_mat")){
  dir.create("../../Intermediate/Count_mat")
}

if(!dir.exists("../../Intermediate/Count_mat/E_coli")){
  dir.create("../../Intermediate/Count_mat/E_coli")
}

# Saving the count matrix and the colData
sampleInfo <- colData(countMatEcoli)
countMatrix <- assay(countMatEcoli)

write.table(sampleInfo, file = "../../Intermediate/Count_mat/E_coli/Sample_data_E_coli.dat")
write.table(countMatrix, file = "../../Intermediate/Count_mat/E_coli/Count_mat_E_coli.dat")
print("Count matrix created for E. coli")

quit(status = 0)
