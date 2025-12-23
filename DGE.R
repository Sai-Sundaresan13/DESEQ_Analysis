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

# Remove unnecessary rows and columns from the counts and the meta file
counts <- counts[ , -c(1:5)]
meta <- meta[ , -c(1:12, 14:28)]

# make sure that the rownames in the metadata file is same as the columnnames in the counts file.
all(rownames(meta) == colnames(counts))

# If the names dont match, change them:
# Method 1: Using Regular Expressions - 
colnames(counts) <- sub("_.*", "", colnames(counts))

# Method 2 - Manually Chnage them
colnames(meta) <- c("GPF", "Stage")

#Method 3 - To assign the rownames of metadata to the column names of the counts or vice verse


#Set one of the attributes as a factor for the DESeq2 analysis so that the DGE is done based on the factor. The factor to be selected is one of the column of the metadata file
#optional before dds creation
meta$Sex <- factor(meta$Sex)   

# Create a DESeq data set object
dds <- DESeqDataSetFromMatrix(countData = counts,
colData = meta,
design = ~ Sex)

dds

# Filter the genes with less counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

# relevel the condition based on the needs
dds$Group <- relevel(dds$Group, ref = "Control")

# Run ESEQ2
dds <- DESeq(dds)
res <- results(dds)
res
summary(res)

# Filter the results based on the p-value
res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)
resultsNames(dds)

# Plot MA plots
plotMA(res)
plotMA(res0.01)

# Order the results and filter it bsased on p-value 
resOrdered=res[order(res$pvalue),]
sourceData=counts(dds,normalized=TRUE)
resFiltered=subset(resOrdered,resOrdered$pvalue<0.01)

# separate the significant genes and plot heatmap
colors=colorRampPalette(c("blue","white","red"))(100)
upRes=resFiltered[resFiltered$log2FoldChange>=2,]
downRes=resFiltered[resFiltered$log2FoldChange<=-2,]

significantUp=sourceData[match((rownames(upRes)),rownames(sourceData)),]
significantDown=sourceData[match((rownames(downRes)),rownames(sourceData)),]

pheatmap(significantUp,scale="row",cluster_rows=FALSE,show_rownames=TRUE,cluster_cols=FALSE,col=colors)
pheatmap(significantDown,scale="row",cluster_rows=FALSE,show_rownames=TRUE,cluster_cols=FALSE,col=colors)

# select the required number of top up and downregulated genes 
topup20=head(upRes,20)
topdown20=head(downRes,20)
top40_combined <- rbind(topup20,topdown20)
significantgenes=sourceData[match((rownames(top40_combined)),rownames(sourceData)),]
pheatmap(significantgenes,scale="row",cluster_rows=FALSE,show_rownames=TRUE,cluster_cols=FALSE,col=colors)



