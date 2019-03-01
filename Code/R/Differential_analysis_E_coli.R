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
library(xtable)

# --------------------------------- Setting up DEseq object analysis -------------------------------------- #
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

# Transforming the data with VTS (n>30), yielding homoskedastic data. TRUE -> unsupervised
# Could also use rld <- rlog(dEcoliFiltered, blind = TRUE), but might be too slow/heavy computationally
dEcoliTransformed <- vst(dEcoliFiltered, blind = TRUE)


# ---------------------------------------------------------------------------------------------------------- #
# ------------------------------------ # Comparing samples # ----------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------- #
colors <- colorRampPalette( rev(brewer.pal(9, "Purples")) )(255)

# Creating the Result/Figure_copy directory if it isn't present. The reason for a second directory is to 
# not overwrite anything in the Result/Figure (since those files are required for Notebook.md file)
if(!dir.exists("../../Results/Figures_copy")){
  dir.create("../../Results/Figures_copy")
}

# ------------------------------------------ # Heat-map # -------------------------------------------------- #
# Function that creates a Poisson distance based Heat-map for all six-samples. The figure is stored in 
# Result/Figures as Pois_dist_heat_E_coli.pdf
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
  
  
  if(exportPdf == TRUE){
    filePath <- "../../Results/Figures/Pois_dist_heat_E_coli.pdf"
    
    # Don't overwrite files used for notebook.md
    if(file.exists(filePath)){
      filePath <- "../../Results/Figures_copy/Pois_dist_heat_E_coli.pdf"
    }
    
    pdf(file = filePath)
  }
  
  if(exportPng == TRUE){
    filePath <- "./../../Results/Figures/Pois_dist_heat_E_coli.png"
    
    # Don't overwrite files used for notebook.md
    if(file.exists(filePath)){
      filePath <- "./../../Results/Figures_copy/Pois_dist_heat_E_coli.png"
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
create_heat_map_pois_dist(dFiltered = dEcoliFiltered, sampleNames = sampleNames, exportPdf = T,exportPng = F)

# ------------------------------------------ # PCA # -------------------------------------------------- #
exportPng = FALSE
exportPdf = TRUE

# Performing a simpel PCA 
pcaData <- plotPCA(dEcoliTransformed, intgroup = c( "condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# For exporting the data
if(exportPdf == TRUE){
  filePath <- "../../Results/Figures/PCA_E_coli.pdf"
  
  # Don't overwrite files used for notebook.md
  if(file.exists(filePath)){
    filePath <- "../../Results/Figures_copy/PCA_E_coli.pdf"
  }
  
  pdf(file = filePath)
}

if(exportPng == TRUE){
  filePath <- "../../Results/Figures/PCA_E_coli.png"
  
  # Don't overwrite files used for notebook.md
  if(file.exists(filePath)){
    filePath <- "../../Results/Figures_copy/PCA_E_coli.png"
  }
  
  png(filename = filePath)
} 

ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) + ggtitle("PCA-plot E. coli") +
  theme(plot.title = element_text(hjust = 0.5)) + geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()

if(exportPdf == TRUE || exportPng == TRUE){
  dev.off()
}


# ---------------------------------------------------------------------------------------------------------- #
# ----------------------------------- # Differential analysis # -------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------- #
# -------------------------------------- # DESeq analysis # ------------------------------------------------ #
DESeqEcoli <- DESeq(dEcoliFiltered)

resultsEcoli <- results(DESeqEcoli, alpha = 0.05)

mcols(resultsEcoli, use.names = T)

summary(resultsEcoli)

# Displaying number of significant p-Values 
nSignPval <- sum(resultsEcoli$padj < 0.05, na.rm = TRUE)
print(sprintf("Number of significant adjusted p-Values = %d", nSignPval))

# Number of upp-regulated and down-regulated significant genes 
nDownReg <- sum((resultsEcoli$padj < 0.05) & (resultsEcoli$log2FoldChange < 0), na.rm = TRUE)
nUpReg <- sum((resultsEcoli$padj < 0.05) & (resultsEcoli$log2FoldChange > 0), na.rm = TRUE)
print( sprintf("Number of upregulated = %d, number of down regulated = %d", nUpReg, nDownReg) )

# ------------------------------ Histogram over p-values ----------------------------------------- #
# Plotting histogram of p-values to see that everything is correct, looks uniform 
# Histogram stored in Results/Figures
exportPng = FALSE
exportPdf = TRUE

# For exporting the data
if(exportPdf == TRUE){
  filePath <- "../../Results/Figures/Histogram_pvalues_E_coli.pdf"
  
  # Don't overwrite files used for notebook.md
  if(file.exists(filePath)){
    filePath <- "../../Results/Figures_copy/Histogram_pvalues_E_coli.pdf"
  }
  
  pdf(file = filePath)
}

if(exportPng == TRUE){
  filePath <- "../../Results/Figures/Histogram_pvalues_E_coli.png"
  
  # Don't overwrite files used for notebook.md
  if(file.exists(filePath)){
    filePath <- "../../Results/Figures_copy/Histogram_pvalues_E_coli.png"
  }
  
  png(filename = filePath)
} 

# Making the plot
hist(resultsEcoli$pvalue[resultsEcoli$baseMean > 1], breaks = 0:20/20,
     col = colors, border = "white")

if(exportPdf == TRUE || exportPng == TRUE){
  dev.off()
}


# ------------------------------------ # Volcano plot # ----------------------------------------------- #
# Add volcano plot
volcPlotGenes <- data.frame(log2FoldChange=resultsEcoli$log2FoldChange, 
                            padj=resultsEcoli$padj)
volcPlotGenes$threshold = as.factor(abs(volcPlotGenes$log2FoldChange) > 2 & volcPlotGenes$padj < 0.05)

# Volcano plot stored in Results/Figures
exportPng = FALSE
exportPdf = TRUE
# For exporting the data

if(exportPdf == TRUE){
  filePath <- "../../Results/Figures/Volcano_E_coli.pdf"
  
  # Don't overwrite files used for notebook.md
  if(file.exists(filePath)){
    filePath <- "../../Results/Figures_copy/Volcano_E_coli.pdf"
  }
  
  pdf(file = filePath)
}

if(exportPng == TRUE){
  filePath <- "../../Results/Figures/Volcano_E_coli.png"
  
  # Don't overwrite files used for notebook.md
  if(file.exists(filePath)){
    filePath <- "../../Results/Figures_copy/Volcano_E_coli.png"
  }
  
  png(filename = filePath)
} 

# Making the plot
ggplot(data=volcPlotGenes, aes(x=log2FoldChange, y=-log10(padj + 10^(-120)), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  theme(legend.position="none") +
  xlim(c(-6, 6)) + ylim(c(0, 125)) +
  xlab("log2 fold change") + ylab("-log10 adj p-value") +
  ggtitle("Volcano plot E.coli") + theme(plot.title = element_text(hjust = 0.5)) 

if(exportPdf == TRUE || exportPng == TRUE){
  dev.off()
}

# ---------------------------------------------------------------------------------------------------------- #
# ------------------------------------- # Exporting results # ---------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------- #
# Exporting result, ordering by p-adj
resultsEcoliOrdered <- resultsEcoli[order(resultsEcoli$padj), ]

# Only exporting significant genes
iToExport <- resultsEcoliOrdered$padj < 0.05
iToExport[is.na(iToExport)] <- FALSE
signResultEcoli <- resultsEcoliOrdered[iToExport, ] 

# Setting correct path for table 
pathToTables <- "../../Results/Tables"
if( !dir.exists(pathToTables) ){
  # The case Tables direcotry doesn't exist
  dir.create(pathToTables)
  pathToExport <- "../../Results/Tables/Table_E_coli.csv"
  rm(pathToTables)
}else if( file.exists("../../Results/Tables/Table_E_coli.csv") ){
  # The case a table alreday exists 
  if( !dir.exists("../../Results/Tables_copy") ){
    dir.create("../../Results/Tables_copy")
  }
  pathToExport <- "../../Results/Tables_copy/Table_E_coli.csv"     
  rm(pathToTables)
}else{
  pathToExport <- "../../Results/Tables/Table_E_coli.csv"
  rm(pathToTables)
}

# Writing result to disk 
write.csv(signResultEcoli, file = pathToExport)

# Creating LaTex table for export, 20 most significant genes 
tableForLatex <- data.frame(log2_fold_change = signResultEcoli$log2FoldChange[1:20], 
                            fold_change = (2 ** (signResultEcoli$log2FoldChange))[1:20],
                            padj = signResultEcoli$padj[1:20], 
                            stringsAsFactors = F)
row.names(tableForLatex) <- row.names(signResultEcoli)[1:20]
codeLatex <- xtable(tableForLatex, digits = -1)
print(codeLatex, file="../../Scratch/LaTex_E_coli_table.txt")
