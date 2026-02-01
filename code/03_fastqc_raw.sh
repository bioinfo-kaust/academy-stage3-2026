#!/bin/bash
# =============================================================================
# Quality Check on Raw Data
# =============================================================================

echo "=============================================="
echo "Step 3: FastQC on Raw Data"
echo "=============================================="

mkdir -p ~/genomics/qc_reports/fastqc_raw
cd ~/genomics

# -----------------------------------------------------------------------------
# Run FastQC on all samples
# -----------------------------------------------------------------------------
echo "Running FastQC on SRR975551..."
fastqc subsampled_data/*.fastq -o qc_reports/fastqc_raw

# -----------------------------------------------------------------------------
# Generate MultiQC summary
# -----------------------------------------------------------------------------
echo "Generating MultiQC report..."
multiqc qc_reports/fastqc_raw -o qc_reports/fastqc_raw -n multiqc_raw
