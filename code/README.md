# RNA-seq Analysis Pipeline

Complete RNA-seq analysis pipeline from raw FASTQ files to biological insights, following the lab instructions from the Bioinformatics Specialization Program.

## Quick Start

```bash
# 1. Setup environment (run once)
./00_setup.sh

# 2. Activate environment
conda activate bioinformatics_course

# 3. Run pipeline with subsampled data (RECOMMENDED for laptops - ~30 min)
./run_pipeline.sh fast

# Or run with full data (slower, several hours)
./run_pipeline.sh all

# Or run individual steps
./run_pipeline.sh download   # Step 1: Download data
./run_pipeline.sh subsample  # Step 1b: Subsample for testing
./run_pipeline.sh qc         # Step 2: Quality check
./run_pipeline.sh trim       # Step 3: Trimming
./run_pipeline.sh align      # Step 4: Alignment
./run_pipeline.sh count      # Step 5: Counting
./run_pipeline.sh de         # Step 6: Differential expression
./run_pipeline.sh viz        # Step 7: Visualization
./run_pipeline.sh enrich     # Step 8: Enrichment
```

## Fast Mode vs Full Mode

| Mode | Command | Reads/Sample | Runtime | Use Case |
|------|---------|--------------|---------|----------|
| Fast | `./run_pipeline.sh fast` | 100,000 | ~30 min | Testing, learning, laptops |
| Full | `./run_pipeline.sh all` with `USE_SUBSAMPLED=false` | All (~20M) | 3-6 hours | Publication, real analysis |

Edit `config.sh` to change `SUBSAMPLE_READS` (default: 100,000) or set `USE_SUBSAMPLED=false` for full data.

## Pipeline Overview

| Step | Script | Description |
|------|--------|-------------|
| 0 | `00_setup.sh` | Create conda environment with all tools |
| 1 | `01_download_data.sh` | Download FASTQ files and reference genome |
| 1b | `01b_subsample_data.sh` | **Subsample reads for faster testing** |
| 2 | `02_quality_check.sh` | Run FastQC, fastp, SeqKit stats |
| 3 | `03_trimming.sh` | Adapter and quality trimming with fastp |
| 4 | `04_alignment.sh` | Align reads with STAR |
| 5 | `05_counting.sh` | Gene quantification with featureCounts |
| 6 | `06_differential_expression.R` | DESeq2 differential expression |
| 7 | `07_visualization.py` | Additional plots with seaborn |
| 8 | `08_enrichment.R` | GO/KEGG enrichment with clusterProfiler |

## Configuration

Edit `config.sh` to customize:
- Sample accessions and conditions
- Number of threads
- Reference genome version
- Analysis parameters

## Default Dataset

The pipeline uses the **Airway dataset (GSE52778)** by default:
- Human airway smooth muscle cells
- Treated with dexamethasone (glucocorticoid)
- 4 samples: 2 control, 2 treated

## Output Structure

```
~/bioinformatics_course/
├── raw_data/           # Downloaded FASTQ files
├── trimmed_data/       # Trimmed reads
├── qc_reports/         # FastQC, fastp, MultiQC reports
├── alignment/          # BAM files and STAR logs
├── counts/             # Count matrices
├── references/         # Genome and annotation files
├── results/
│   ├── tables/         # DE results, normalized counts
│   ├── figures/        # All visualizations
│   └── enrichment/     # GO/KEGG results
└── logs/               # Pipeline logs
```

## Requirements

- 8GB RAM minimum (16GB recommended)
- 20GB disk space
- Internet connection (for data download)
- Unix-like environment (Linux, macOS, WSL)

## Troubleshooting

**Conda not found:**
```bash
source ~/.bashrc  # or ~/.zshrc
```

**Permission denied:**
```bash
chmod +x *.sh *.R *.py
```

**STAR index fails (memory):**
- Use `GENOME_SUBSET="chr22"` in config.sh for testing
- Need ~32GB RAM for full human genome

**R packages fail:**
```R
BiocManager::install("DESeq2", force = TRUE)
```
