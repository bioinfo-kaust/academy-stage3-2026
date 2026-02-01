#!/bin/bash
# =============================================================================
# Install Conda Environment
# =============================================================================

echo "=============================================="
echo "Setting up Bioinformatics Environment"
echo "=============================================="

# Create conda environment
mamba create -n genomics python=3.11 -y
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate genomics

# Install core bioinformatics tools using mamba (faster than conda)
mamba install -c conda-forge -c bioconda sra-tools fastqc fastp seqkit multiqc star salmon samtools igv -y

# Install Python data science packages
mamba install -c conda-forge -c bioconda pandas numpy matplotlib seaborn scikit-learn openpyxl gprofiler-official -y

# Test that tools are accessible
fastqc --version
fastp --version
STAR --version
samtools --version
salmon --version
multiqc --version

# Test Python packages
python -c "import pandas; import seaborn; print('Python packages OK')"


# Install R and its packages using mamba

mamba install -c conda-forge -c bioconda \
    r-base=4.4 \
    r-tidyverse \
    r-biocmanager \
    bioconductor-deseq2 \
    bioconductor-tximport \
    bioconductor-clusterprofiler \
    bioconductor-org.hs.eg.db \
    bioconductor-enrichplot \
    r-pheatmap \
    r-ggrepel \
    r-ggplot2 \
	r-pheatmap \
	r-rcolorbrewer

# Test R packages
R -q -e "suppressWarnings({library(ggplot2); library(pheatmap); library(RColorBrewer); library(DESeq2); library(tximport); library(clusterProfiler)}); cat('All R packages loaded successfully\n')"


echo "=============================================="
echo "Setup complete!"
echo "Activate with: conda activate genomics2"
echo "=============================================="
