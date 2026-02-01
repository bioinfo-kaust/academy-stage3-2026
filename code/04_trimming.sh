#!/bin/bash
# =============================================================================
# Trim Adapters and Low Quality Bases
# =============================================================================

echo "=============================================="
echo "Step 4: Trimming with fastp"
echo "=============================================="

mkdir -p ~/genomics/trimmed_data
mkdir -p ~/genomics/qc_reports/fastp
cd ~/genomics

# -----------------------------------------------------------------------------
# Trim All samples - for loop
# -----------------------------------------------------------------------------
for SAMPLE_ID in SRR975551 SRR975552 SRR975553 SRR975554 SRR975555 SRR975556; do
    echo "Trimming ${SAMPLE_ID}..."
    fastp \
        --in1 subsampled_data/${SAMPLE_ID}_1.fastq \
        --in2 subsampled_data/${SAMPLE_ID}_2.fastq \
        --out1 trimmed_data/${SAMPLE_ID}_1.trimmed.fastq \
        --out2 trimmed_data/${SAMPLE_ID}_2.trimmed.fastq \
        --qualified_quality_phred 20 \
        --length_required 36 \
        --detect_adapter_for_pe \
        --thread 4 \
        --json qc_reports/fastp/${SAMPLE_ID}.json \
        --html qc_reports/fastp/${SAMPLE_ID}.html
done

# -----------------------------------------------------------------------------
# Generate MultiQC summary
# -----------------------------------------------------------------------------
echo "Generating MultiQC report..."
multiqc qc_reports/fastp -o qc_reports/fastp -n multiqc_fastp
