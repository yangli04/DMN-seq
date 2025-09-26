[![DOI](https://zenodo.org/badge/965719200.svg)](https://doi.org/10.5281/zenodo.17211256)

# DMN-seq Data Processing Pipeline

This repository contains comprehensive data processing workflows for DMN-seq analysis. The pipeline provides two distinct analytical approaches for methylation profiling:

1. **5mC Detection** (`5mC_detection/`) - Quantitative detection and analysis of 5-methylcytosine modifications
2. **Hypomethylation Enrichment** (`hypomethylation_enrichment/`) - Identification and characterization of hypomethylated regions

## Installation

### Prerequisites
- Python 3.10
- Snakemake 7.32.4
- Mamba or Conda package manager (Mamba recommended for faster dependency resolution)

### Environment Setup

Create and activate the required environment:

```bash
# Using mamba (recommended)
mamba env create -f environment.yml
mamba activate dmnseq

# Or using conda
conda env create -f environment.yml
conda activate dmnseq
```

## Configuration

Before running the analysis, configure the pipeline parameters:

1. **Update workflow configuration**: Edit `config.yaml` in your chosen analysis directory
2. **Configure reference genome paths**: Specify genome assembly and annotation files
3. **Set software paths**: Update paths to required bioinformatics tools
4. **Prepare input data**: Place raw sequencing files in the `raw_data/` directory

Example configuration files are provided in each analysis subdirectory.

## Usage

### Option 1: Cluster Execution (Recommended for HPC)

The pipeline is optimized for SLURM-based high-performance computing environments (tested on University of Chicago Midway3 cluster):

```bash
cd hypomethylation_enrichment/  # or 5mC_detection/
snakemake --profile ../snakemake_config
```

**Note**: Modify `./snakemake_config/config.yaml` to match your cluster specifications.

### Option 2: Local Execution

For local analysis or smaller datasets:

```bash
cd hypomethylation_enrichment/  # or 5mC_detection/
snakemake --cores 12  # adjust core count as appropriate
```

## Workflow Overview

### Hypomethylation Enrichment Analysis
Identifies and analyzes regions with reduced methylation levels, suitable for:
- Hypomethylated region detection
- Differential methylation analysis
- Methylation landscape characterization

![Hypomethylation Enrichment Workflow](hypomethylation_enrichment/dag.svg)

### 5mC Detection Analysis
Quantitative analysis of 5-methylcytosine modifications, optimized for:
- Single-base resolution methylation calling
- Quantitative methylation profiling
- Context-specific methylation analysis

![5mC Detection Workflow](5mC_detection/dag.svg)

## Output Structure

The pipeline generates structured output directories containing:
- Quality control reports
- Processed alignment files
- Methylation call files
- Statistical analysis results
- Visualization plots and summaries

## Citation

- cite this software

  ```BibTex
  @software{yangli_2025_dmnseq_zenodo.17211256,
    author       = {Yang Li},
    title        = {v1.0.0},
    publisher = {Zenodo},
    year         = {2025},
    doi          = {10.5281/zenodo.17211256},
    url          = {https://doi.org/10.5281/zenodo.17211256},
  }

  ```

## Support

For questions, issues, and comments, please refer to the repository's issue tracker or contact __Yang Li__ (U Chicago) yliuchicago@uchicago.edu. 
