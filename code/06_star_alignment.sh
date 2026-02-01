#!/bin/bash
# =============================================================================
# Genome Alignment with STAR (OPTIONAL)
# =============================================================================

echo "=============================================="
echo "Step 6: STAR Alignment (Optional)"
echo "=============================================="

mkdir -p ~/genomics/references/star_index
mkdir -p ~/genomics/alignment
cd ~/genomics

# -----------------------------------------------------------------------------
# Build STAR index
# -----------------------------------------------------------------------------
echo "Building STAR index... - once per species or genome reference file set"
STAR --runMode genomeGenerate \
    --runThreadN 4 \
    --genomeDir references/star_index \
    --genomeFastaFiles references/Homo_sapiens.GRCh38.dna.chromosome.22.fa \
    --sjdbGTFfile references/Homo_sapiens.GRCh38.110.chr22.gtf \
    --sjdbOverhang 100 \
    --genomeSAindexNbases 11

# -----------------------------------------------------------------------------
# Perform Alignment for All Samples
# -----------------------------------------------------------------------------

for SAMPLE_ID in SRR975551 SRR975552 SRR975553 SRR975554 SRR975555 SRR975556; do
    echo "Aligning ${SAMPLE_ID}..."
    STAR --runThreadN 4 \
        --genomeDir references/star_index \
        --readFilesIn trimmed_data/${SAMPLE_ID}_1.trimmed.fastq trimmed_data/${SAMPLE_ID}_2.trimmed.fastq \
        --outFileNamePrefix alignment/${SAMPLE_ID}_ \
        --outSAMtype BAM SortedByCoordinate
done

# -----------------------------------------------------------------------------
# Index BAM files
# -----------------------------------------------------------------------------
samtools index alignment/*.bam

