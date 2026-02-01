# =============================================================================
# 09_deseq2.R - Differential Expression Analysis
# =============================================================================
# Run with: Rscript 09_deseq2.R
# =============================================================================

library(DESeq2)
library(tximport)

cat("==============================================\n")
cat("Step 9: DESeq2 Differential Expression\n")
cat("==============================================\n")

# -----------------------------------------------------------------------------
# Create output directories
# -----------------------------------------------------------------------------
dir.create("~/genomics/results", showWarnings = FALSE)
dir.create("~/genomics/results/tables", showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
cat("Loading data...\n")

txi <- readRDS("~/genomics/counts/tximport.rds")
sample_info <- read.table("~/genomics/counts/sample_info.tsv", header = TRUE)
sample_info$condition <- factor(sample_info$condition, levels = c("Normal", "Tumor"))

cat(paste0("  Samples: ", nrow(sample_info), "\n"))
cat(paste0("  Genes: ", nrow(txi$counts), "\n"))

# -----------------------------------------------------------------------------
# Create DESeq2 dataset
# -----------------------------------------------------------------------------
cat("Running DESeq2...\n")

dds <- DESeqDataSetFromTximport(txi, colData = sample_info, design = ~ condition)

# Filter low count genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
cat(paste0("  Genes after filtering: ", nrow(dds), "\n"))

# Run DESeq2
dds <- DESeq(dds)

# Get results (Tumor vs Normal)
results <- results(dds, contrast = c("condition", "Tumor", "Normal"))
results <- results[order(results$padj), ]

# -----------------------------------------------------------------------------
# Filter significant genes
# -----------------------------------------------------------------------------
sig_genes <- subset(results, padj < 0.05 & abs(log2FoldChange) > 1)

cat(paste0("\nSignificant genes (padj<0.05, |log2FC|>1): ", nrow(sig_genes), "\n"))
cat(paste0("  Upregulated in Tumor: ", sum(sig_genes$log2FoldChange > 0, na.rm = TRUE), "\n"))
cat(paste0("  Downregulated in Tumor: ", sum(sig_genes$log2FoldChange < 0, na.rm = TRUE), "\n"))

# -----------------------------------------------------------------------------
# Save results
# -----------------------------------------------------------------------------
cat("\nSaving results...\n")

# All results
all_results <- as.data.frame(results)
all_results$gene_id <- rownames(all_results)
all_results <- all_results[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
write.csv(all_results, "~/genomics/results/tables/deseq2_all_results.csv", row.names = FALSE)

# Significant only
sig_results <- as.data.frame(sig_genes)
sig_results$gene_id <- rownames(sig_results)
sig_results <- sig_results[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
write.csv(sig_results, "~/genomics/results/tables/deseq2_significant.csv", row.names = FALSE)

# Normalized counts
norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
norm_counts$gene_id <- rownames(norm_counts)
write.csv(norm_counts, "~/genomics/results/tables/normalized_counts.csv", row.names = FALSE)

# DESeq2 object
saveRDS(dds, "~/genomics/results/deseq2_object.rds")

cat("==============================================\n")
cat("DESeq2 complete!\n")
cat("Outputs:\n")
cat("  - results/tables/deseq2_all_results.csv\n")
cat("  - results/tables/deseq2_significant.csv\n")
cat("  - results/tables/normalized_counts.csv\n")
cat("==============================================\n")
