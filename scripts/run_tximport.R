#!/usr/bin/env Rscript
# =============================================================================
# run_tximport.R - Aggregate Salmon Transcript Counts to Gene Level
# =============================================================================
# This script:
#   1. Generates transcript-to-gene mapping from a GTF file
#   2. Imports Salmon quantification results
#   3. Aggregates transcript counts to gene level
#   4. Outputs count matrices for DESeq2
#
# Usage:
#   Rscript run_tximport.R --gtf <gtf_file> --salmon_dir <salmon_quant_dir> --outdir <output_dir>
#
# Example:
#   Rscript run_tximport.R \
#       --gtf ~/genomics/references/Homo_sapiens.GRCh38.110.chr11.gtf \
#       --salmon_dir ~/genomics/salmon_quant \
#       --outdir ~/genomics/counts
#
# Required R packages: tximport, optparse
# Install with: install.packages(c("tximport", "optparse"))
# =============================================================================

# Load required packages
suppressPackageStartupMessages({
    library(optparse)
    library(tximport)
})

# -----------------------------------------------------------------------------
# Parse command-line arguments
# -----------------------------------------------------------------------------
option_list <- list(
    make_option(c("-g", "--gtf"), type = "character", default = NULL,
                help = "Path to GTF annotation file [required]",
                metavar = "FILE"),
    make_option(c("-s", "--salmon_dir"), type = "character", default = NULL,
                help = "Path to Salmon quantification directory containing sample subdirectories [required]",
                metavar = "DIR"),
    make_option(c("-o", "--outdir"), type = "character", default = "counts",
                help = "Output directory for count matrices [default: %default]",
                metavar = "DIR"),
    make_option(c("-t", "--tx2gene"), type = "character", default = NULL,
                help = "Path to existing tx2gene.tsv file (optional, will be generated from GTF if not provided)",
                metavar = "FILE"),
    make_option(c("--samples"), type = "character", default = NULL,
                help = "Comma-separated list of sample names (optional, auto-detected from salmon_dir if not provided)",
                metavar = "SAMPLES"),
    make_option(c("-v", "--verbose"), action = "store_true", default = TRUE,
                help = "Print progress messages [default: %default]")
)

opt_parser <- OptionParser(
    option_list = option_list,
    usage = "Usage: %prog [options]",
    description = "\nAggregate Salmon transcript counts to gene level using tximport.\nGenerates tx2gene mapping from GTF if not provided."
)

opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$salmon_dir)) {
    stop("Error: --salmon_dir is required. Use --help for usage information.")
}

if (is.null(opt$gtf) && is.null(opt$tx2gene)) {
    stop("Error: Either --gtf or --tx2gene must be provided. Use --help for usage information.")
}

# -----------------------------------------------------------------------------
# Helper function: Generate tx2gene from GTF
# -----------------------------------------------------------------------------
generate_tx2gene_from_gtf <- function(gtf_file, verbose = TRUE) {
    if (verbose) cat("Generating tx2gene mapping from GTF...\n")
    if (verbose) cat(paste0("  Reading: ", gtf_file, "\n"))

    # Read GTF file
    gtf <- read.table(gtf_file, sep = "\t", header = FALSE, comment.char = "#",
                      stringsAsFactors = FALSE, quote = "")

    # Filter for transcript features
    gtf_transcripts <- gtf[gtf$V3 == "transcript", ]

    if (nrow(gtf_transcripts) == 0) {
        stop("Error: No transcript features found in GTF file.")
    }

    if (verbose) cat(paste0("  Found ", nrow(gtf_transcripts), " transcripts\n"))

    # Extract transcript_id and gene_id from attributes (column 9)
    extract_attribute <- function(attributes, attr_name) {
        pattern <- paste0(attr_name, ' "([^"]+)"')
        matches <- regmatches(attributes, regexec(pattern, attributes))
        sapply(matches, function(m) if (length(m) > 1) m[2] else NA)
    }

    tx2gene <- data.frame(
        transcript_id = extract_attribute(gtf_transcripts$V9, "transcript_id"),
        gene_id = extract_attribute(gtf_transcripts$V9, "gene_id"),
        stringsAsFactors = FALSE
    )

    # Remove any NAs
    tx2gene <- tx2gene[complete.cases(tx2gene), ]

    # Remove version numbers from transcript IDs (e.g., ENST00000123456.1 -> ENST00000123456)
    tx2gene$transcript_id <- gsub("\\.\\d+$", "", tx2gene$transcript_id)

    if (verbose) cat(paste0("  Generated mapping for ", nrow(tx2gene), " transcripts\n"))

    return(tx2gene)
}

# -----------------------------------------------------------------------------
# Main workflow
# -----------------------------------------------------------------------------
if (opt$verbose) {
    cat("==============================================\n")
    cat("tximport: Aggregate Transcript Counts to Gene Level\n")
    cat("==============================================\n\n")
}

# Step 1: Get or generate tx2gene mapping
if (!is.null(opt$tx2gene) && file.exists(opt$tx2gene)) {
    if (opt$verbose) cat("Loading existing tx2gene mapping...\n")
    tx2gene <- read.table(opt$tx2gene, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    colnames(tx2gene) <- c("transcript_id", "gene_id")
    if (opt$verbose) cat(paste0("  Loaded ", nrow(tx2gene), " transcript-gene pairs\n"))
} else {
    if (is.null(opt$gtf)) {
        stop("Error: GTF file required to generate tx2gene mapping.")
    }
    if (!file.exists(opt$gtf)) {
        stop(paste0("Error: GTF file not found: ", opt$gtf))
    }
    tx2gene <- generate_tx2gene_from_gtf(opt$gtf, opt$verbose)
}

# Step 2: Find sample directories in salmon_dir
if (opt$verbose) cat("\nFinding Salmon quantification files...\n")

if (!is.null(opt$samples)) {
    # Use provided sample list
    samples <- unlist(strsplit(opt$samples, ","))
} else {
    # Auto-detect samples from directory
    sample_dirs <- list.dirs(opt$salmon_dir, full.names = FALSE, recursive = FALSE)
    # Filter to only directories containing quant.sf
    samples <- sample_dirs[sapply(sample_dirs, function(s) {
        file.exists(file.path(opt$salmon_dir, s, "quant.sf"))
    })]
}

if (length(samples) == 0) {
    stop("Error: No Salmon quantification files found in ", opt$salmon_dir)
}

if (opt$verbose) cat(paste0("  Found ", length(samples), " samples: ", paste(samples, collapse = ", "), "\n"))

# Build file paths
files <- file.path(opt$salmon_dir, samples, "quant.sf")
names(files) <- samples

# Verify all files exist
missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
    stop("Error: Missing quant.sf files:\n  ", paste(missing_files, collapse = "\n  "))
}

# Step 3: Import with tximport
if (opt$verbose) cat("\nImporting Salmon counts with tximport...\n")

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

if (opt$verbose) cat(paste0("  Imported ", nrow(txi$counts), " genes across ", ncol(txi$counts), " samples\n"))

# Step 4: Create output directory
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# Step 5: Save outputs
if (opt$verbose) cat("\nSaving output files...\n")

# Gene count matrix (rounded for DESeq2)
counts <- as.data.frame(round(txi$counts))
counts$gene_id <- rownames(counts)
counts <- counts[, c("gene_id", samples)]
counts_file <- file.path(opt$outdir, "gene_counts.tsv")
write.table(counts, counts_file, sep = "\t", quote = FALSE, row.names = FALSE)
if (opt$verbose) cat(paste0("  ", counts_file, "\n"))

# TPM matrix (for visualization)
tpm <- as.data.frame(txi$abundance)
tpm$gene_id <- rownames(tpm)
tpm <- tpm[, c("gene_id", samples)]
tpm_file <- file.path(opt$outdir, "gene_tpm.tsv")
write.table(tpm, tpm_file, sep = "\t", quote = FALSE, row.names = FALSE)
if (opt$verbose) cat(paste0("  ", tpm_file, "\n"))

# Sample info (auto-generate based on sample names)
# Assumes sample names start with condition (e.g., KO_1, WT_1)
conditions <- gsub("_[0-9]+$", "", samples)
sample_info <- data.frame(
    sample_id = samples,
    condition = conditions,
    stringsAsFactors = FALSE
)
sample_info_file <- file.path(opt$outdir, "sample_info.tsv")
write.table(sample_info, sample_info_file, sep = "\t", quote = FALSE, row.names = FALSE)
if (opt$verbose) cat(paste0("  ", sample_info_file, "\n"))

# Save tx2gene mapping (for reference)
tx2gene_file <- file.path(opt$outdir, "tx2gene.tsv")
write.table(tx2gene, tx2gene_file, sep = "\t", quote = FALSE, row.names = FALSE)
if (opt$verbose) cat(paste0("  ", tx2gene_file, "\n"))

# Save tximport object for DESeq2
txi_file <- file.path(opt$outdir, "tximport.rds")
saveRDS(txi, txi_file)
if (opt$verbose) cat(paste0("  ", txi_file, "\n"))

# Print summary
if (opt$verbose) {
    cat("\n==============================================\n")
    cat("Summary\n")
    cat("==============================================\n")
    cat(paste0("Samples processed: ", length(samples), "\n"))
    cat(paste0("Genes quantified:  ", nrow(txi$counts), "\n"))
    cat(paste0("Output directory:  ", opt$outdir, "\n"))
    cat("\nOutput files:\n")
    cat("  gene_counts.tsv  - Raw counts for DESeq2\n")
    cat("  gene_tpm.tsv     - TPM values for visualization\n")
    cat("  sample_info.tsv  - Sample metadata\n")
    cat("\nDone!\n")
}
