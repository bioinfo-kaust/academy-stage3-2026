#!/bin/bash
# =============================================================================
# run_all.sh - Run Complete RNA-seq Pipeline
# =============================================================================
# Dataset: GSE50760 - Colorectal Cancer (Normal vs Tumor)
# =============================================================================

set -e
cd "$(dirname "$0")"

echo "=============================================="
echo "Running Complete RNA-seq Pipeline"
echo "=============================================="

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

echo "=============================================="
echo "Pipeline complete!"
echo "Results: ~/rnaseq_course/results/"
echo "=============================================="
