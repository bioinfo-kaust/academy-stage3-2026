# =============================================================================
# Aggregate Transcript Counts to Gene Level
# =============================================================================
# Run with: Rscript 08_tximport.R
# =============================================================================

library(tximport)

# -----------------------------------------------------------------------------
# Create output directory
# -----------------------------------------------------------------------------
dir.create("~/genomics/counts", showWarnings = FALSE)

cat("Loading tx2gene mapping from file...\n")
tx2gene <- read.table("~/genomics/references/transcriptID2geneID.tsv",
                  sep = "\t", header = TRUE, comment.char = "")
colnames(tx2gene) <- c("transcript_id", "gene_id")

cat(paste0("  Loaded ", nrow(tx2gene), " transcript-gene pairs\n"))
# -----------------------------------------------------------------------------
# Define samples
# -----------------------------------------------------------------------------
samples <- c("SRR975551", "SRR975552", "SRR975553",
             "SRR975554", "SRR975555", "SRR975556")
conditions <- c("Normal", "Normal", "Normal", "Tumor", "Tumor", "Tumor")

# -----------------------------------------------------------------------------
# Import Salmon counts
# -----------------------------------------------------------------------------
cat("Importing Salmon counts...\n")

quant_files <- c(
  "~/genomics/salmon_quant/SRR975551/quant.sf",
  "~/genomics/salmon_quant/SRR975552/quant.sf",
  "~/genomics/salmon_quant/SRR975553/quant.sf",
  "~/genomics/salmon_quant/SRR975554/quant.sf",
  "~/genomics/salmon_quant/SRR975555/quant.sf",
  "~/genomics/salmon_quant/SRR975556/quant.sf"
)
names(quant_files) <- samples

txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

cat(paste0("  Imported ", nrow(txi$counts), " genes\n"))

# -----------------------------------------------------------------------------
# Save outputs
# -----------------------------------------------------------------------------
cat("Saving outputs...\n")

# Count matrix
counts <- as.data.frame(round(txi$counts))
counts$gene_id <- rownames(counts)
counts <- counts[, c("gene_id", samples)]
write.table(counts, "~/genomics/counts/count_matrix.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# TPM matrix
tpm <- as.data.frame(txi$abundance)
tpm$gene_id <- rownames(tpm)
tpm <- tpm[, c("gene_id", samples)]
write.table(tpm, "~/genomics/counts/tpm_matrix.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Sample info
sample_info <- data.frame(sample_id = samples, condition = conditions)
write.table(sample_info, "~/genomics/counts/sample_info.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# tximport object for DESeq2
saveRDS(txi, "~/genomics/counts/tximport.rds")

