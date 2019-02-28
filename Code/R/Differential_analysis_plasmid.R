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

# ---------------------------------------------------------------------------------------------------------- #
# ----------------------------------- # Differential analysis # -------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------- #
# -------------------------------------- # DESeq analysis # ------------------------------------------------ #
DESeqPlasmid <- DESeq(dPlasmidFiltered)

resultsPlasmid <- results(DESeqPlasmid, alpha = 0.05)

summary(resultsPlasmid)

# Displaying number of significant p-Values 
nSignPval <- sum(resultsPlasmid$padj < 0.05, na.rm = TRUE)
print(sprintf("Number of significant adjusted p-Values = %d", nSignPval))

# Number of upp-regulated and down-regulated significant genes 
nDownReg <- sum((resultsPlasmid$padj < 0.05) & (resultsPlasmid$log2FoldChange < 0), na.rm = TRUE)
nUpReg <- sum((resultsPlasmid$padj < 0.05) & (resultsPlasmid$log2FoldChange > 0), na.rm = TRUE)
print( sprintf("Number of upregulated = %d, number of down regulated = %d", nUpReg, nDownReg) )

# Plotting the top-gene
topGene <- rownames(resultsPlasmid)[which.min(resultsPlasmid$padj)]
plotCounts(DESeqPlasmid, gene = topGene, intgroup=c("condition"))

# ------------------------------ Histogram over p-values ----------------------------------------- #
# Histogram stored in Results/Figures
exportPng = F
exportPdf = F
# For exporting the data

if(exportPdf == TRUE){
  filePath <- "../../Results/Figures/Histogram_pvalues_plasmid.pdf"
  
  # Don't overwrite files used for notebook.md
  if(file.exists(filePath)){
    filePath <- "../../Results/Figures_copy/Histogram_pvalues_plasmid.pdf"
  }
  
  pdf(file = filePath)
}

if(exportPng == TRUE){
  filePath <- "../../Results/Figures/Histogram_pvalues_plasmid.png"
  
  # Don't overwrite files used for notebook.md
  if(file.exists(filePath)){
    filePath <- "../../Results/Figures_copy/Histogram_pvalues_plasmid.png"
  }
  
  png(filename = filePath)
} 

hist(resultsPlasmid$pvalue[resultsPlasmid$baseMean > 1], breaks = 0:20/20,
     col = "steelblue", border = "white", main = "P-values plasmid")



if(exportPdf == TRUE || exportPng == TRUE){
  dev.off()
}

rm(exportPdf, exportPng, filePath)

# ------------------------------------ # Volcano plot # ----------------------------------------------- #
# Volcano plot stored in Results/Figures
exportPng = F; exportPdf = F

if(exportPdf == TRUE){
  filePath <- "../../Results/Figures/Volcano_Plasmid.pdf"
  
  # Don't overwrite files used for notebook.md
  if(file.exists(filePath)){
    filePath <- "../../Results/Figures_copy/Volcano_Plasmid.pdf"
  }
  
  pdf(file = filePath)
}

if(exportPng == TRUE){
  filePath <- "../../Results/Figures/Volcano_Plasmid.png"
  
  # Don't overwrite files used for notebook.md
  if(file.exists(filePath)){
    filePath <- "../../Results/Figures_copy/Volcano_Plasmid.png"
  }
  
  png(filename = filePath)
} 

volcPlotGenes <- data.frame(log2FoldChange=resultsPlasmid$log2FoldChange, 
                            padj=resultsPlasmid$padj)
volcPlotGenes$threshold = as.factor(abs(volcPlotGenes$log2FoldChange) > 0.5 & volcPlotGenes$padj < 0.05)

# Making the plot
ggplot(data=volcPlotGenes, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  theme(legend.position="none") +
  xlim(c(-2.5, 2.5)) + ylim(c(0, 15)) +
  xlab("log2 fold change") + ylab("-log10 adj p-value") +
  ggtitle("Volcano-plot plasmid") + theme(plot.title = element_text(hjust = 0.5)) 


if(exportPdf == TRUE || exportPng == TRUE){
  dev.off()
}

rm(exportPdf, exportPng, volcPlotGenes, filePath)
# ------------------------------------ # Exporting result # ------------------------------------------ #
# Function that given a string will find the corresponding match in geneInfo$Gene_names
# Input:
# strInput, string to search for in geneInfo$Gene_name (geneInfoName)
# Output
# index if matching, else it returns -1
find_index_match <- function(strInput, geneInfoName)
{
  for(a in 1:length(geneInfoName)){
    if(geneInfoName[a] == strInput){
      return(a)
    }
  }
  
  return(-1) 
}

# Function thatt will annotate the table with descriptions
# Input:
# geneInfo from the file Annotation_plasmid.tsv in the data table
# Data-frame with result from DESeq
# Ouput:
# Annoated table (gene descriptions)
annotate_result <- function(geneInfo, tableToExport)
{
  for( i in 1:dim(tableToExport)[1] ){
    geneName <- rownames(tableToExport)[i]
    indexGene <- find_index_match(geneName, geneInfo$Gene_name)
    #    print(indexGene)
    if( indexGene  == -1){
      next
    }else{
      tableToExport$Description[i] <- geneInfo$Description[indexGene]
    }
  }
  
  return(tableToExport)
}

# Only export sorted result
resultsPlasmidOrdered <- resultsPlasmid[order(resultsPlasmid$padj), ]

# Read the gene names and description 
pathToAnnotation = "../../Data/Reference_data/Plasmid/Annotation_plasmid.tsv"
if(!file.exists(pathToAnnotation)){
  print("Annotation file doesn't exist")
  quit(status = 1)
}

geneInfo <- read.table(pathToAnnotation, header = T, sep = "\t")
geneInfo <- data.frame(lapply(geneInfo, as.character), stringsAsFactors=FALSE)

# Data frame to export 
tableToExport <- data.frame(baseMean=resultsPlasmid$baseMean, log2FoldChange=resultsPlasmid$log2FoldChange, 
                            lfcSE=resultsPlasmid$lfcSE, stat=resultsPlasmid$stat, 
                            padj=resultsPlasmid$padj, Description=character(length(resultsPlasmid$baseMean)), 
                            stringsAsFactors = F)
rownames(tableToExport) <- rownames(resultsPlasmid)

# Annotate the table 
tableToExport <- annotate_result(geneInfo, tableToExport)






