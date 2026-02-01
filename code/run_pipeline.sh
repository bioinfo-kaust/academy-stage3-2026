#!/bin/bash
# =============================================================================
# RNA-seq Analysis Pipeline Overview
# =============================================================================
# Dataset: GSE50760 - Colorectal Cancer
# Normal: SRR975551, SRR975552, SRR975553
# Tumor:  SRR975554, SRR975555, SRR975556
# =============================================================================

cat << 'EOF'
==============================================
RNA-seq Analysis Pipeline
==============================================

Pipeline Steps:

  00_setup.sh          - Install conda environment (run once)
  01_download_data.sh  - Download data from SRA
  02_subsample_data.sh - Subsample to 1M reads
  03_fastqc_raw.sh     - Quality check (raw data)
  04_trimming.sh       - Trim adapters and low quality
  05_fastqc_trimmed.sh - Quality check (trimmed data)
  06_star_alignment.sh - STAR alignment (OPTIONAL)
  07_salmon_quant.sh   - Salmon quantification
  08_tximport.R        - Gene-level counts
  09_deseq2.R          - Differential expression
  10_visualization.py  - Create plots
  11_enrichment.R      - GO/KEGG enrichment

How to run:

  # First time setup:
  conda activate base
  bash 00_setup.sh

  # Run pipeline:
  conda activate rnaseq_course
  bash 01_download_data.sh
  bash 02_subsample_data.sh
  bash 03_fastqc_raw.sh
  bash 04_trimming.sh
  bash 05_fastqc_trimmed.sh
  bash 07_salmon_quant.sh
  Rscript 08_tximport.R
  Rscript 09_deseq2.R
  python3 10_visualization.py
  Rscript 11_enrichment.R

  # Or run everything at once:
  bash run_all.sh

==============================================
EOF
