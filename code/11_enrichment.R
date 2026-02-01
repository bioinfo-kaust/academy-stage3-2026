# =============================================================================
# 11_enrichment.R - GO and KEGG Enrichment Analysis
# =============================================================================
# Run with: Rscript 11_enrichment.R
# =============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

cat("==============================================\n")
cat("Step 11: Enrichment Analysis\n")
cat("==============================================\n")

# -----------------------------------------------------------------------------
# Create output directory
# -----------------------------------------------------------------------------
dir.create("~/genomics/results/enrichment", showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Load significant genes
# -----------------------------------------------------------------------------
cat("Loading significant genes...\n")

sig_genes <- read.csv("~/genomics/results/tables/deseq2_significant.csv")
cat(paste0("  Found ", nrow(sig_genes), " significant genes\n"))

# Remove version numbers from gene IDs
all_genes <- gsub("\\..*", "", sig_genes$gene_id)
up_genes <- gsub("\\..*", "", sig_genes$gene_id[sig_genes$log2FoldChange > 0])
down_genes <- gsub("\\..*", "", sig_genes$gene_id[sig_genes$log2FoldChange < 0])

cat(paste0("  Upregulated: ", length(up_genes), "\n"))
cat(paste0("  Downregulated: ", length(down_genes), "\n"))

# -----------------------------------------------------------------------------
# Convert Ensembl to Entrez IDs
# -----------------------------------------------------------------------------
cat("Converting gene IDs...\n")

gene_mapping <- bitr(all_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat(paste0("  Mapped ", nrow(gene_mapping), " genes\n"))

up_entrez <- gene_mapping$ENTREZID[gene_mapping$ENSEMBL %in% up_genes]
down_entrez <- gene_mapping$ENTREZID[gene_mapping$ENSEMBL %in% down_genes]

# -----------------------------------------------------------------------------
# GO Enrichment - Upregulated genes
# -----------------------------------------------------------------------------
cat("\nRunning GO enrichment...\n")

if (length(up_entrez) >= 3) {
  go_up <- enrichGO(gene = up_entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                    pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

  if (!is.null(go_up) && nrow(go_up) > 0) {
    cat(paste0("  GO terms (upregulated): ", nrow(go_up), "\n"))
    write.csv(as.data.frame(go_up), "~/genomics/results/enrichment/go_upregulated.csv", row.names = FALSE)

    pdf("~/genomics/results/enrichment/go_upregulated.pdf", width = 10, height = 8)
    print(dotplot(go_up, showCategory = 15, title = "GO: Upregulated in Tumor"))
    dev.off()
  }
}

# -----------------------------------------------------------------------------
# GO Enrichment - Downregulated genes
# -----------------------------------------------------------------------------
if (length(down_entrez) >= 3) {
  go_down <- enrichGO(gene = down_entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                      pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

  if (!is.null(go_down) && nrow(go_down) > 0) {
    cat(paste0("  GO terms (downregulated): ", nrow(go_down), "\n"))
    write.csv(as.data.frame(go_down), "~/genomics/results/enrichment/go_downregulated.csv", row.names = FALSE)

    pdf("~/genomics/results/enrichment/go_downregulated.pdf", width = 10, height = 8)
    print(dotplot(go_down, showCategory = 15, title = "GO: Downregulated in Tumor"))
    dev.off()
  }
}

# -----------------------------------------------------------------------------
# KEGG Pathway Analysis
# -----------------------------------------------------------------------------
cat("Running KEGG enrichment...\n")

all_entrez <- gene_mapping$ENTREZID

if (length(all_entrez) >= 3) {
  kegg <- enrichKEGG(gene = all_entrez, organism = "hsa",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05)

  if (!is.null(kegg) && nrow(kegg) > 0) {
    cat(paste0("  KEGG pathways: ", nrow(kegg), "\n"))
    write.csv(as.data.frame(kegg), "~/genomics/results/enrichment/kegg_pathways.csv", row.names = FALSE)

    pdf("~/genomics/results/enrichment/kegg_pathways.pdf", width = 10, height = 8)
    print(dotplot(kegg, showCategory = 15, title = "KEGG Pathways"))
    dev.off()
  }
}

cat("==============================================\n")
cat("Enrichment analysis complete!\n")
cat("Results in: ~/genomics/results/enrichment/\n")
cat("==============================================\n")
