#!/bin/bash
# =============================================================================
# Transcript Quantification with Salmon
# =============================================================================

echo "=============================================="
echo "Step 7: Salmon Quantification"
echo "=============================================="

mkdir -p ~/genomics/salmon_index
mkdir -p ~/genomics/salmon_quant
cd ~/genomics

# -----------------------------------------------------------------------------
# Build Salmon index
# -----------------------------------------------------------------------------
echo "Building Salmon index... - only once per species and genome reference"
salmon index \
    -t references/transcriptome_chr22.fa \
    -i references/salmon_index \
    --threads 4


# -----------------------------------------------------------------------------
# Perform Salmon Quant for All Samples
# -----------------------------------------------------------------------------

for SAMPLE_ID in SRR975551 SRR975552 SRR975553 SRR975554 SRR975555 SRR975556; do
    echo "Quantifying ${SAMPLE_ID}..."
    salmon quant \
        -i references/salmon_index \
        -l A \
        -1 trimmed_data/${SAMPLE_ID}_1.trimmed.fastq \
        -2 trimmed_data/${SAMPLE_ID}_2.trimmed.fastq \
        -o salmon_quant/${SAMPLE_ID} \
        --threads 4 \
        --validateMappings \
        --gcBias \
        --seqBias
done
