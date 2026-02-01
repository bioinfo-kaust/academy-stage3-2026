#!/bin/bash
# =============================================================================
# Install Conda Environment
# =============================================================================

echo "=============================================="
echo "Setting up Bioinformatics Environment"
echo "=============================================="

# Create conda environment
mamba create -n genomics2 python=3.9 -y
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate genomics2

# Install bioinformatics tools
echo "Installing bioinformatics tools..."
mamba install -c conda-forge -c bioconda sra-tools fastqc fastp seqkit multiqc star salmon samtools -y

# Install Python packages
echo "Installing Python packages..."
pip install pandas numpy matplotlib seaborn scikit-learn

# Install R and packages
echo "Installing R..."
mamba install \
  -c conda-forge -c bioconda \
  r-base=4.4 \
  r-tidyverse \
  r-biocmanager \
  bioconductor-deseq2 \
  bioconductor-tximport \
  bioconductor-clusterprofiler \
  bioconductor-org.hs.eg.db \
  bioconductor-enrichplot \
  r-pheatmap \
  r-ggrepel

echo "=============================================="
echo "Setup complete!"
echo "Activate with: conda activate genomics2"
echo "=============================================="
