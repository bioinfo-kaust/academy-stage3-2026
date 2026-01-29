# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an **educational curriculum project** for a 4-day Bioinformatics Specialization Program (Feb 4-7, 2026) focused on RNA-sequencing analysis. The goal is to create a static GitHub Pages website presenting the course content.

## Project Structure

- `program-info.txt` - Complete curriculum: learning outcomes, 4-day schedule, session descriptions, labs
- `instructions.txt` - Website requirements and specifications
- `slides/day[1-4].pptx` - Placeholder presentation files (to be replaced with actual content)

## Website Deliverables

The task is to create a static website with:
1. **Home page** - Program overview, learning outcomes, practical info
2. **Schedule page** - Daily timetables for all 4 days with session/topic details
3. **Labs index page** - List of all labs linking to individual lab pages
4. **Individual lab pages** - Detailed instructions with bash/Python/R code examples
5. **Environment setup page** - Linux terminal and Conda setup instructions (including WSL for Windows)

## Technical Context

### Bioinformatics Tools Referenced
- **QC/Preprocessing:** FastQC, fastp, TrimGalore, Cutadapt
- **Alignment:** STAR, samtools, IGV
- **Quantification:** Salmon, StringTie, featureCounts
- **Differential Expression:** DESeq2 (R)
- **Visualization:** seaborn (Python), MultiQC
- **Enrichment Analysis:** ClusterProfiler (R)

### Environment Requirements for Lab Content
- Linux terminal (WSL for Windows users)
- Conda with bioconda channel
- Python 3.9+ with seaborn
- R with DESeq2 and ClusterProfiler

## Key Guidelines

- Target audience: Students with little computational background
- Use step-by-step instructions with command examples
- Use code blocks for all shell/Python/R commands
- Keep language student-friendly and instructional
- Labs should use consistent datasets throughout (preprocessing → analysis → enrichment)
- Website must be responsive and have clear navigation
