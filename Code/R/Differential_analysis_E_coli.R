# This file performs the differential analysis of for the E.coli genome 
# using the DEseq2-package. Output from this file is:
# Tables (can be found in Result/tables): 
# Table with most significant genes 
# Figures (can be found in Result/Figures):
# Heat-map of difference between samples 
rm(list = ls())

library(DESeq2)
library(GenomicAlignments)
library(Rsamtools)
library(GenomicFeatures)

# --------------------------------- Differential analysis -------------------------------------- #
# File path to count matrix and sample data
pathCountMat <- "../../Intermediate/Count_mat/E_coli/Count_mat_E_coli.dat"
pathSampleData <- "../../Intermediate/Count_mat/E_coli/Sample_data_E_coli.dat"

# If data doesn't exist exit program 
if(!file.exists(pathCountMat) || !file.exists(pathSampleData)){
  message("The count matrix or sample data doesn't exist for E.coli")
  quit(status = 1)
}

# Reading the data
countMatEcoli <- read.table(file = pathCountMat)
sampleDataEcoli <- read.table(file = pathSampleData)

# Setting dose values from NA to 0
sampleDataEcoli$dose[1:3] <- c(0L, 0L, 0L)

# Setting none imipenem as base-line
sampleDataEcoli$condition <- relevel(sampleDataEcoli$condition, "none")

# Creating DEseq2 object 
dEcoli <- DESeqDataSetFromMatrix(countData = countMatEcoli, 
                                 colData = sampleDataEcoli, 
                                 design = ~ condition)

# Pre-filtering the data set by removing rows with no counts,
# or only a single count across all samples. (9 rows)
dEcoliFiltered <- dEcoli[ rowSums(counts(dEcoli)) > 1, ]

# Transforming the data with VTS (n>30), yielding homoskedastic data
# Could also use rld <- rlog(dEcoliFiltered, blind = FALSE), but might be too slow/heavy computationally
dEcoliTransformed <- vst(dEcoliFiltered, blind = FALSE)

# Creating heatmap of sample-to-sample distances
sampleDists <- dist(t(assay(dEcoliTransformed)))

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( dEcoliTransformed$condition, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)


