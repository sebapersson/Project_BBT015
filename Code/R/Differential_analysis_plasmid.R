# This file performs the differential analysis of for the plasmid
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
library(pheatmap)
library(PoiClaClu)
library(ggplot2)
library(RColorBrewer)

# --------------------------------- Setting up DEseq object analysis -------------------------------------- #
# File path to count matrix and sample data
pathCountMat <- "../../Intermediate/Count_mat/Plasmid/Count_mat_plasmid.dat"
pathSampleData <- "../../Intermediate/Count_mat/Plasmid/Sample_data_plasmid.dat"

# If data doesn't exist exit program 
if(!file.exists(pathCountMat) || !file.exists(pathSampleData)){
  message("The count matrix or sample data doesn't exist for E.coli")
  quit(status = 1)
}

# Reading the data
countMatPlasmid <- read.table(file = pathCountMat)
sampleDataPlasmid <- read.table(file = pathSampleData)

# Setting dose values from NA to 0
sampleDataPlasmid$dose[1:3] <- c(0L, 0L, 0L)

# Setting none imipenem as base-line
sampleDataPlasmid$condition <- relevel(sampleDataEcoli$condition, "none")

# Creating DEseq2 object 
dPlasmid <- DESeqDataSetFromMatrix(countData = countMatPlasmid, 
                                   colData = sampleDataPlasmid, 
                                   design = ~ condition)

# Pre-filtering the data set by removing rows with no counts,
# or only a single count across all samples. (2 rows)
dPlasmidFiltered <- dPlasmid[ rowSums(counts(dPlasmid)) > 1, ]

# Normalising the data (used for PCA with other things)
# blind = TRUE -> unsupervised transformation
dPlasmidTransformed <- varianceStabilizingTransformation(dPlasmidFiltered, blind = TRUE)
# ---------------------------------------------------------------------------------------------------------- #
# ------------------------------------ # Comparing samples # ----------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------- #

# Creating the Result/Figure_copy directory if it isn't present. The reason for a second directory is to 
# not overwrite anything in the Result/Figure (since those files are required for Notebook.md file)
if(!dir.exists("../../Results/Figures_copy")){
  dir.create("../../Results/Figures_copy")
}

# ------------------------------------------ # Heat-map # -------------------------------------------------- #
# Function that creates a Poisson distance based Heat-map for all six-samples. The figure is stored in 
# Result/Figures as Pois_dist_heat_plasmid.pdf
# Input:
# DEseq filtered object
# Sample names for the heat-map 
# Output
# Figure 
create_heat_map_pois_dist <- function(dFiltered, sampleNames, exportPdf=FALSE, exportPng=FALSE)
{
  # Calculate the Poisson distance
  poisDist <- PoissonDistance(t(counts(dFiltered)))
  samplePoisDistMatrix <- as.matrix( poisDist$dd )
  rownames(samplePoisDistMatrix) <- sampleNames
  colnames(samplePoisDistMatrix) <- NULL
  
  # Choose colour purple 
  colors <- colorRampPalette( rev(brewer.pal(9, "Purples")) )(255)
  
  # For exporting the data
  if(exportPdf == TRUE){
    filePath <- "../../Results/Figures/Pois_dist_heat_plasmid.pdf"
    
    # Don't overwrite files used for notebook.md
    if(file.exists(filePath)){
      filePath <- "../../Results/Figures_copy/Pois_dist_heat_plasmid.pdf"
    }
    
    pdf(file = filePath)
  }
  
  if(exportPng == TRUE){
    filePath <- "./../../Results/Figures/Pois_dist_heat_plasmid.png"
    
    # Don't overwrite files used for notebook.md
    if(file.exists(filePath)){
      filePath <- "./../../Results/Figures_copy/Pois_dist_heat_plasmid.png"
    }
    
    png(filename = filePath)
  }
  
  pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisDist$dd,
         clustering_distance_cols = poisDist$dd,
         col = colors)
  
  if(exportPdf == TRUE || exportPng == TRUE){
    dev.off()
  }

}

# Creat the Poisson heat-map
sampleNames <- c("Sample1-Cont.", "Sample2-Cont.", "Sample3-Cont.", 
                 "Sample4-Case", "Sample5-Case" ,"Sample6-Case")
create_heat_map_pois_dist(dFiltered = dPlasmidFiltered, sampleNames = sampleNames, exportPdf  =  F)

# ------------------------------------------ # Heat-map # -------------------------------------------------- #
# This section was formally a function, however using a function seemed to clash with the plot-creation, 
# Hence this is now a continous code-section

# Chose what format to export 
exportPng = F; exportPdf = F

# Performing a simpel PCA 
pcaData <- plotPCA(dPlasmidTransformed, intgroup = c( "condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# For exporting the data
if(exportPdf == TRUE){
  filePath <- "../../Results/Figures/PCA_plasmid.pdf"
  
  # Don't overwrite files used for notebook.md
  if(file.exists(filePath)){
    filePath <- "../../Results/Figures_copy/PCA_plasmid.pdf"
  }
  
  pdf(file = filePath)
}

if(exportPng == TRUE){
  filePath <- "./../../Results/Figures/PCA_plasmid.png"
  
  # Don't overwrite files used for notebook.md
  if(file.exists(filePath)){
    filePath <- "./../../Results/Figures_copy/PCA_plasmid.png"
  }
  
  png(filename = filePath)
} 

ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) + ggtitle("PCA-plot plasmid") +
  theme(plot.title = element_text(hjust = 0.5)) + geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()

if(exportPdf == TRUE || exportPng == TRUE){
  dev.off()
}
  
# Remove uneccesary variables 
rm(exportPdf, exportPng, percentVar, pcaData, filePath)






