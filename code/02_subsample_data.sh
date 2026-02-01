#!/bin/bash
# =============================================================================
# Subsample to 1M reads per sample
# =============================================================================
# Makes pipeline faster for teaching purposes
# =============================================================================

echo "=============================================="
echo "Step 2: Subsample Data (1M reads)"
echo "=============================================="

mkdir -p ~/genomics/subsampled_data
cd ~/genomics

# -----------------------------------------------------------------------------
# Subsample for All Files - run in a loop
# -----------------------------------------------------------------------------
for SAMPLE_ID in SRR975551 SRR975552 SRR975553 SRR975554 SRR975555 SRR975556; do
    echo "Subsampling ${SAMPLE_ID}..."
    fastp \
        --in1 raw_data/${SAMPLE_ID}_1.fastq \
        --in2 raw_data/${SAMPLE_ID}_2.fastq \
        --out1 subsampled_data/${SAMPLE_ID}_1.fastq \
        --out2 subsampled_data/${SAMPLE_ID}_2.fastq \
        --reads_to_process 1000000 \
        --disable_quality_filtering \
        --disable_length_filtering \
        --disable_adapter_trimming \
        --json /dev/null --html /dev/null
done
