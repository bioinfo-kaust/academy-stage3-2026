#!/bin/bash
# =============================================================================
# Quality Check on Trimmed Data
# =============================================================================

echo "=============================================="
echo "Step 5: FastQC on Trimmed Data"
echo "=============================================="

mkdir -p ~/genomics/qc_reports/fastqc_trimmed
cd ~/genomics

# -----------------------------------------------------------------------------
# Run FastQC on all trimmed samples
# -----------------------------------------------------------------------------
echo "Running FastQC on all files"
fastqc trimmed_data/*.fastq -o qc_reports/fastqc_trimmed

# -----------------------------------------------------------------------------
# Generate MultiQC summary
# -----------------------------------------------------------------------------
echo "Generating MultiQC report..."
multiqc qc_reports/fastqc_trimmed -o qc_reports/fastqc_trimmed -n multiqc_trimmed
