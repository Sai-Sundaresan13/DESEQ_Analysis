###############################################################################
# DESeq2 Analysis Pipeline - R Script
#
# Consolidates:
#   - DESeq2 differential expression analysis (based on original DGE.R)
#   - Visualization: Volcano plot, MA plot, PCA plot, Heatmap (per README.md)
#
# Input: the count matrix produced by pipeline.sh (featureCounts output),
#        plus a sample metadata (colData) file you provide.
###############################################################################

# ============================== CONFIG ======================================
counts_file      <- "counts/counts.txt"   # featureCounts output (tab-delimited)
metadata_file    <- "meta.csv"            # sample metadata, rownames = sample IDs
condition_column <- "Group"               # column in metadata to test on
reference_level  <- "Control"             # baseline/control level for that column
padj_cutoff      <- 0.05
lfc_cutoff       <- 1                     # |log2FoldChange| threshold for "significant"
top_n_heatmap    <- 20                    # top up/down genes per side in final heatmap
output_dir       <- "results"
# =============================================================================

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 0. Load required packages ---------------------------------------------
suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
    library(pheatmap)
    library(dplyr)
})

# ---- 1. Load counts and metadata --------------------------------------------
# featureCounts output has 6 leading annotation columns:
# Geneid, Chr, Start, End, Strand, Length -- keep Geneid as rownames, drop the rest.
raw_counts <- read.table(counts_file, header = TRUE, row.names = 1,
                          sep = "\t", comment.char = "#", check.names = FALSE)
counts <- raw_counts[, -(1:5)]

# Clean up featureCounts' long BAM-path column names, e.g.
# "aligned/sample1_sorted.bam" -> "sample1"
colnames(counts) <- gsub(".*/", "", colnames(counts))
colnames(counts) <- gsub("_sorted\\.bam$", "", colnames(counts))

meta <- read.csv(metadata_file, header = TRUE, row.names = 1)

# ---- 2. Match samples between counts and metadata ---------------------------
common <- intersect(colnames(counts), rownames(meta))
if (length(common) == 0) {
    stop("No overlapping sample names between counts and metadata. ",
         "Check colnames(counts) vs rownames(meta).")
}
counts <- counts[, common]
meta   <- meta[common, , drop = FALSE]
stopifnot(all(colnames(counts) == rownames(meta)))

meta[[condition_column]] <- factor(meta[[condition_column]])

# ---- 3. Build DESeq2 dataset -------------------------------------------------
design_formula <- as.formula(paste("~", condition_column))
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData   = meta,
                               design    = design_formula)

# Filter low-count genes
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]

# Set the reference/control level for the comparison
dds[[condition_column]] <- relevel(dds[[condition_column]], ref = reference_level)

# ---- 4. Run DESeq2 -----------------------------------------------------------
dds <- DESeq(dds)
res <- results(dds, alpha = padj_cutoff)

cat("\n=== DESeq2 summary ===\n")
summary(res)

resOrdered <- res[order(res$pvalue), ]
write.csv(as.data.frame(resOrdered), file.path(output_dir, "DESeq2_results.csv"))

# Normalized counts, used by all the plots below
normCounts <- counts(dds, normalized = TRUE)

# ---- 5. Significant genes ----------------------------------------------------
resDF <- as.data.frame(res)
resDF$gene <- rownames(resDF)
resDF <- resDF[!is.na(resDF$padj), ]

sigUp   <- resDF %>% filter(padj < padj_cutoff, log2FoldChange >=  lfc_cutoff)
sigDown <- resDF %>% filter(padj < padj_cutoff, log2FoldChange <= -lfc_cutoff)

write.csv(sigUp,   file.path(output_dir, "upregulated_genes.csv"),   row.names = FALSE)
write.csv(sigDown, file.path(output_dir, "downregulated_genes.csv"), row.names = FALSE)

cat(sprintf("\nSignificant genes: %d up, %d down (padj < %s, |log2FC| >= %s)\n",
            nrow(sigUp), nrow(sigDown), padj_cutoff, lfc_cutoff))

# ============================ VISUALIZATION ==================================

# ---- a. Volcano Plot ---------------------------------------------------------
resDF$sig <- "NS"
resDF$sig[resDF$padj < padj_cutoff & resDF$log2FoldChange >=  lfc_cutoff] <- "Up"
resDF$sig[resDF$padj < padj_cutoff & resDF$log2FoldChange <= -lfc_cutoff] <- "Down"

volcano <- ggplot(resDF, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c(Up = "red", Down = "blue", NS = "grey70")) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    theme_minimal() +
    labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10(padj)", color = "")

ggsave(file.path(output_dir, "volcano_plot.png"), volcano, width = 7, height = 6, dpi = 300)

# ---- b. MA Plot ---------------------------------------------------------------
png(file.path(output_dir, "MA_plot.png"), width = 1400, height = 1200, res = 200)
plotMA(res, alpha = padj_cutoff, main = "MA Plot")
dev.off()

# ---- c. PCA Plot ----------------------------------------------------------
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = condition_column, returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = .data[[condition_column]])) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_minimal() +
    labs(title = "PCA Plot")

ggsave(file.path(output_dir, "pca_plot.png"), pca_plot, width = 7, height = 6, dpi = 300)

# ---- d. Heatmap of top up/down regulated genes ------------------------------
colors <- colorRampPalette(c("blue", "white", "red"))(100)

topUp   <- head(sigUp[order(sigUp$padj), "gene"],   top_n_heatmap)
topDown <- head(sigDown[order(sigDown$padj), "gene"], top_n_heatmap)
topGenes <- c(topUp, topDown)

if (length(topGenes) > 0) {
    heatmapData <- normCounts[topGenes, , drop = FALSE]

    png(file.path(output_dir, "heatmap_top_genes.png"), width = 1600, height = 1800, res = 200)
    pheatmap(heatmapData, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
             show_rownames = TRUE, col = colors,
             main = paste("Top", length(topUp), "up /", length(topDown), "down regulated genes"))
    dev.off()
} else {
    message("No significant genes found for heatmap at current thresholds.")
}

cat("\nAll results and plots written to: ", normalizePath(output_dir), "\n")
