#!/bin/bash
# =============================================================================
# Download Data from SRA
# =============================================================================
# Dataset: GSE50760 - Colorectal Cancer RNA-seq
# Normal: SRR975551, SRR975552, SRR975553
# Tumor:  SRR975554, SRR975555, SRR975556
# =============================================================================

echo "=============================================="
echo "Step 1: Download Data from SRA"
echo "=============================================="

# Create directories
mkdir -p ~/genomics/raw_data
mkdir -p ~/genomics/references
cd ~/genomics/raw_data

# -----------------------------------------------------------------------------
# Download Normal sample 1 or see next command to download all in parallel
# -----------------------------------------------------------------------------
prefetch --progress SRR975551

# -----------------------------------------------------------------------------
# Download All Files in Parallel - to make it faster [Note: Requires large disk space ~ 30 GB]
# -----------------------------------------------------------------------------
printf "%s\n" SRR975551 SRR975552 SRR975553 SRR975554 SRR975555 SRR975556 \
| xargs -n 1 -P 6 prefetch --progress
fasterq-dump --split-files --threads 8 --progress SRR*/*.sra
# -----------------------------------------------------------------------------
# Convert the downloaded sra files into fastq using fasterq-dump [Note: Requires large disk space ~ 140 GB]
# -----------------------------------------------------------------------------
fasterq-dump --split-files --threads 8 --progress SRR*/*.sra

# -----------------------------------------------------------------------------
# To save time and disk space, you can downloaded already downsampled datasets below
# -----------------------------------------------------------------------------
link to google drive gzip file containing all subset_fastq files [subsampled_data]

# -----------------------------------------------------------------------------
# Download reference genome (chromosome 22 only)
# -----------------------------------------------------------------------------
echo "Downloading reference files..."
cd ~/genomics/references

wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
gunzip Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz

wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
gunzip Homo_sapiens.GRCh38.110.gtf.gz
grep -E "^#|^22	" Homo_sapiens.GRCh38.110.gtf > Homo_sapiens.GRCh38.110.chr22.gtf

wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip -c Homo_sapiens.GRCh38.cdna.all.fa.gz | awk '/^>/ {keep = /chromosome:GRCh38:22:/} keep' > transcriptome_chr22.fa

#Generate Transcript - Gene ID mapping file by parsing the GTF file
awk -F'\t' '$3=="transcript" {
  if (match($9, /gene_id "[^"]+"/)) {
    gene = substr($9, RSTART, RLENGTH)
    sub(/^gene_id "/, "", gene)
    sub(/"$/, "", gene)
  }
  if (match($9, /transcript_id "[^"]+"/)) {
    tx = substr($9, RSTART, RLENGTH)
    sub(/^transcript_id "/, "", tx)
    sub(/"$/, "", tx)
  }
  if (gene && tx) print tx "\t" gene}' Homo_sapiens.GRCh38.110.chr22.gtf > transcriptID2geneID.tsv
# -----------------------------------------------------------------------------