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
library(HDF5Array)

# --------------------------------- Differential analysis -------------------------------------- #
# File path to summarised experiments
pathCountMat <- "../../Intermediate/Count_mat/E_coli/"

# If data doesn't exist exit program 
if(!dir.exists(pathCountMat)){
  message("The count matrix doesn't exist for E.coli")
  quit(status = 1)
}

# Reading the data
countMatEcoli <- loadHDF5SummarizedExperiment(pathCountMat)

assay(countMatEcoli)
