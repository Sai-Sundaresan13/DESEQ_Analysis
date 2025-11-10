#This file contains the script for Differential gene Expression analysis using DESeq2 Pakage in R.

# Load the required Packages
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(dplyr)
library(corrplot)


# load the counts file and the metadata file
counts <- read.csv("counts1.csv", header = TRUE, row.names = 1)
meta <- read.csv("meta1.csv", header = TRUE, row.names = 1)

# make sure that the rownames in the metadata file is same as the columnnames in the counts file.
all(rownames(meta) == colnames(counts))

#Set one of the attributes as a factor for the DESeq2 analysis so that the DGE is done based on the factor. The factor to be selected is one of the column of the metadata file
meta$Sex <- factor(meta$Sex)   

# Create a DESeq data set object
dds <- DESeqDataSetFromMatrix(countData = counts,
colData = meta,
design = ~ Sex)

