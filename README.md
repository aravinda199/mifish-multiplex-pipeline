# MiFish eDNA Multiplex Pipeline

An R-based bioinformatics pipeline for processing environmental DNA (eDNA) metabarcoding data generated using multiplexed MiFish primers (MiFish-U and MiFish-E). 

This workflow utilizes `dada2` for denoising and ASV inference, and `cutadapt` for primer removal.

## Pipeline Overview
1. Quality profiling and ambiguous base (N) filtering.
2. Primer detection and trimming using `cutadapt` (handles both MiFish-U and MiFish-E in forward/reverse orientations).
3. Quality filtering and trimming (`truncLen=c(110, 105)`).
4. Error rate learning and ASV inference via `dada2` (pseudo-pooling).
5. Merging of paired-end reads and chimera removal.
6. Output generation: ASV FASTA, ASV abundance table, and read-tracking statistics.

## Prerequisites

### 1. R Packages
You will need the following libraries installed in your R environment:
* `dada2`
* `ShortRead`
* `Biostrings`
* `ggplot2`

### 2. Cutadapt
This pipeline requires [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) to be installed on your system. 

## Usage Instructions

1. Clone this repository to your local machine.
2. Place your raw paired-end sequencing files (`*_R1.fastq.gz` and `*_R2.fastq.gz`) into the `data/raw/` directory.
3. Open `scripts/mifish_dada2_pipeline.R`.
4. **Important:** Update the `CUTADAPT_PATH` variable at the top of the script to match the installation path of cutadapt on your system.
5. Run the script. Output files will be generated in the `output/` directory.

Maintained by: Ara | Last Updated: 09.09.2026
