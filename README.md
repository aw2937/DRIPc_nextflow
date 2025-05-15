# DRIP-Seq Single-Read Processing Pipeline

This Nextflow pipeline processes single-read DRIP-seq (DNA-RNA Immunoprecipitation sequencing) data. It automates steps from raw FASTQ files to alignment, coverage generation, peak calling, and annotation.

## Table of Contents

1.  [Overview](#overview)
2.  [Features](#features)
3.  [Requirements](#requirements)
4.  [Setup](#setup)
    *   [Pipeline Code](#pipeline-code)
    *   [Input Data](#input-data)
    *   [Reference Genomes & Indices](#reference-genomes--indices)
    *   [Custom Scripts & Models](#custom-scripts--models)
5.  [Usage](#usage)
    *   [Input Samplesheet](#input-samplesheet)
    *   [Running the Pipeline](#running-the-pipeline)
6.  [Parameters](#parameters)
    *   [Mandatory Parameters](#mandatory-parameters)
    *   [Optional Parameters](#optional-parameters)
    *   [Reference Genome Parameters (hg19 default)](#reference-genome-parameters-hg19-default)
    *   [Custom Script Paths](#custom-script-paths)
7.  [Output Directory Structure](#output-directory-structure)
8.  [Dependency Management](#dependency-management)
9. [Troubleshooting](#troubleshooting)

## Overview

The pipeline takes raw single-end FASTQ files, performs quality control, trimming, alignment, and various downstream analyses typical for DRIP-seq, including:
*   Merging per-lane FASTQ files for each sample.
*   Quality control using FastQC and MultiQC.
*   Adapter and quality trimming using Trim Galore.
*   Alignment to a reference genome (hg19 by default) using Bowtie2.
*   BAM file processing: sorting, duplicate marking (using samtools fixmate and markdup), and indexing.
*   Generation of strand-specific coverage files (bedGraph and BigWig).
*   Splitting BAM files into forward and reverse strands.
*   Optional: Peak calling using MACS3.
*   Optional: Peak annotation using a custom DROPA script.
*   Optional: Alternative peak calling using the "Chedin method" (samtools depth, custom perl scripts, stochHMM).

## Features

*   **Automated Workflow:** End-to-end processing from FASTQ to downstream analysis.
*   **Parallelization:** Leverages Nextflow's parallel execution capabilities for efficient processing of multiple samples.
*   **Reproducibility:** Ensures consistent analysis by defining specific software versions (via Conda) and parameters.
*   **Modularity:** Optional steps like MACS3, DROPA, and Chedin method can be skipped.
*   **Standardized Outputs:** Generates common bioinformatics file formats.

## Requirements

*   **Nextflow:** (version 21.10.x or later recommended).
*   **Conda:** To manage software dependencies. Alternatively, Docker or Singularity can be configured if you adapt the process definitions.
*   **Required Tools:** The pipeline uses Conda to install the following (among others):
    *   `fastqc`
    *   `multiqc`
    *   `trim-galore` (and its dependency `cutadapt`)
    *   `bowtie2`
    *   `samtools`
    *   `bedtools`
    *   `ucsc-wigtobigwig`
    *   `macs3` (optional)
    *   `python3` (with pandas, numpy, scipy for DROPA - optional)
    *   `perl` (for Chedin method custom scripts - optional)
    *   `stochhmm` (for Chedin method - optional, **Note:** this might require manual installation or a specific container as it's not commonly available in standard Conda channels).

## Setup

### 1. Pipeline Code
Clone or download this repository/pipeline script (`dripseq.nf`).

```bash
# If it's in a git repository:
# git clone <repository_url>
# cd <repository_directory>


### 2. Input Data
*   Place your gzipped FASTQ files in a directory. The path to this directory will be specified by `params.input_dir`.
*   The pipeline uses a samplesheet (see [Input Samplesheet](#input-samplesheet)) to find these files.

### 3. Reference Genomes & Indices
You need to provide paths to:
*   **Bowtie2 genome index prefix:** e.g., `/path/to/ref_genome/hg19/hg19` (where files like `hg19.1.bt2`, `hg19.2.bt2` etc. exist).
*   **Genome FASTA index file (`.fai`):** e.g., `/path/to/ref_genome/hg19/hg19.fa.fai`.
    These are specified via parameters like `params.bowtie2_index_hg19` and `params.genome_fai_hg19`.

### 4. Custom Scripts & Models
Certain steps use custom scripts or models. By default, the pipeline expects these to be in a `bin/` directory relative to the `dripseq.nf` script. You can override these paths using parameters.
*   `bin/DROPA_v1.0.0.py`: For peak annotation.
*   `bin/normalize.pl`: For Chedin method WIG normalization.
*   `bin/wig2fa.pl`: For Chedin method WIG to FASTA conversion.
*   `bin/DRIPc.hmm`: HMM model for stochHMM in the Chedin method.

Ensure the Perl scripts (`normalize.pl`, `wig2fa.pl`) are executable:
```bash
chmod +x bin/normalize.pl bin/wig2fa.pl
```

Additionally, for DROPA, you'll need:
*   DROPA gene reference directory (e.g., `GeneReference/hg19_RefSeq/`) specified by `params.dropa_ref_hg19`.
*   DROPA genome size file (e.g., `GeneReference/hg19.genome`) specified by `params.dropa_gsize_hg19`.

## Usage

### 1. Input Samplesheet
Create a CSV file (e.g., `samples.csv`) with the following header and structure:

```csv
sample,fastq_glob
EH-1H,EH-1H_S1_L00*_R1_001.fastq.gz
EH-2H,EH-2H_S2_L00*_R1_001.fastq.gz
EH-3H,EH-3H_S3_L00*_R1_001.fastq.gz
# ... and so on for all your samples
```
*   `sample`: A unique identifier for your sample. This will be used for naming output files.
*   `fastq_glob`: A glob pattern relative to `params.input_dir` that matches the FASTQ file(s) for that sample. The pipeline will `cat` all matching files for a sample.

### 2. Running the Pipeline
Execute the pipeline using the `nextflow run` command:

```bash
nextflow run dripseq.nf \
    --samplesheet samples.csv \
    --input_dir /path/to/your/fastq/directory \
    --outdir ./results \
    --bowtie2_index_hg19 /path/to/ref_genome/hg19/hg19 \
    --genome_fai_hg19 /path/to/ref_genome/hg19/hg19.fa.fai \
    # Add other parameters as needed, e.g., to enable optional steps
    # --skip_macs3 false
    # --skip_dropa false
    # --skip_chedin false
```

You can also use a `nextflow.config` file to specify parameters and process configurations (e.g., executor, CPU, memory).

## Parameters

### Mandatory Parameters
*   `--samplesheet`: Path to the input samplesheet CSV file (e.g., `samples.csv`).
*   `--input_dir`: Path to the directory containing raw FASTQ files (e.g., `data/fastq`).
*   `--bowtie2_index_hg19`: Path to the Bowtie2 index prefix for hg19 (e.g., `ref_genome/hg19/hg19`).
*   `--genome_fai_hg19`: Path to the genome FASTA index file for hg19 (e.g., `ref_genome/hg19/hg19.fa.fai`).

### Optional Parameters
*   `--outdir`: Directory where results will be published (Default: `results`).
*   `--skip_fastqc`: Skip FastQC steps (Default: `false`).
*   `--skip_multiqc`: Skip MultiQC report generation (Default: `false`).
*   `--skip_macs3`: Skip MACS3 peak calling (Default: `false`).
*   `--skip_dropa`: Skip DROPA peak annotation (Default: `true`).
*   `--skip_chedin`: Skip the Chedin peak calling method (Default: `true`).

### Reference Genome Parameters (hg19 default)
*   `--bowtie2_index_hg19`: Bowtie2 index prefix for hg19.
*   `--genome_fai_hg19`: Genome FASTA index for hg19.
*   `--macs_gsize_hg19`: MACS3 genome size for hg19 (Default: `hs`).
*   `--dropa_ref_hg19`: Path to DROPA gene reference directory for hg19 (Default: `GeneReference/hg19_RefSeq/`).
*   `--dropa_gsize_hg19`: Path to DROPA genome size file for hg19 (Default: `GeneReference/hg19.genome`).

*Note: If analyzing genomes other than hg19, you would need to adapt these parameters and potentially the workflow logic if tools behave differently.*

### Custom Script Paths
These default to locations within a `bin/` directory relative to the pipeline script.
*   `--dropa_script`: Path to `DROPA_v1.0.0.py` (Default: `${projectDir}/bin/DROPA_v1.0.0.py`).
*   `--normalize_script`: Path to `normalize.pl` (Default: `${projectDir}/bin/normalize.pl`).
*   `--wig2fa_script`: Path to `wig2fa.pl` (Default: `${projectDir}/bin/wig2fa.pl`).
*   `--dripc_hmm_model`: Path to `DRIPc.hmm` (Default: `${projectDir}/bin/DRIPc.hmm`).

## Output Directory Structure

The pipeline will create an output directory (specified by `params.outdir`, default `results/`) with the following structure:

```
results/
├── 01_merged_fastq/        # Merged FASTQ files per sample
├── 02_fastqc_raw/          # FastQC reports for merged FASTQs
├── 03_trimmed_fastq/       # Trimmed FASTQ files and TrimGalore reports
├── 04_fastqc_trimmed/      # FastQC reports for trimmed FASTQs
├── 05_bam_files/           # Processed BAM files (aligned, sorted, duplicates marked) and BAI indices
├── 06_coverage_bigwig/     # Strand-specific BigWig coverage files
├── 07_split_bams/          # BAM files split by strand (fwd/rev) with BAI indices
├── 08_macs3_peaks/         # (If not skipped) MACS3 peak calling results per sample
│   └── <sample_id>/
├── 09_dropa_annotation/    # (If not skipped) DROPA peak annotation results per sample
│   └── <sample_id>/
├── 10_chedin_method/       # (If not skipped) Results from the Chedin peak calling method
│   ├── custom_fa/
│   ├── normalized_wigs/
│   ├── stochhmm_peaks/
│   └── wigs/
└── multiqc/                # MultiQC report aggregating QC results
    ├── multiqc_report.html
    └── multiqc_data/
```

## Dependency Management

Software dependencies are managed using Conda. Each process in the Nextflow script specifies its required Conda environment. Nextflow will automatically create (or use cached) these environments.
Ensure Conda is installed and accessible in your `PATH`.

For tools not easily available via Conda (like `stochhmm`), you might need to:
1.  Install it manually into a Conda environment that Nextflow can use.
2.  Modify the process to use a Docker/Singularity container where the tool is pre-installed.
3.  Ensure the tool is available in your system `PATH` and remove the `conda` directive for that specific process (less reproducible).


## Troubleshooting

*   **Custom Script Not Found:** Ensure `DROPA_v1.0.0.py`, `normalize.pl`, `wig2fa.pl`, `DRIPc.hmm` are in the correct location (default: `bin/` directory relative to `dripseq.nf`) or that you've specified the correct paths using parameters (`--dropa_script`, etc.). Make sure Perl scripts are executable.
*   **`stochhmm` not found:** This tool can be difficult to install. Consider using a container or ensuring it's in a Conda environment accessible to Nextflow.
*   **`wigToBigWig` errors:** This tool sometimes has issues with Conda installation. Ensure the `ucsc-wigtobigwig` package is correctly installed.
*   **Conda environment creation issues:** Check your Conda setup and internet connection. You can try creating the environments manually first.
*   **File/Path Issues:** Double-check all input file paths, especially for reference genomes and indices. Nextflow is case-sensitive and requires exact paths.
```
