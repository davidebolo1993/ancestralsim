# ancestralsim

Simulate ancient DNA (aDNA) reads from pangenome haplotypes

## Overview

`ancestralsim` is a pipeline for simulating ancient DNA sequencing reads from locus-specific pangenome haplotypes. It uses [gargammel](https://github.com/grenaud/gargammel) to generate realistic aDNA reads with authentic damage patterns, fragmentation, and contamination profiles, then aligns them to a reference chromosome using BWA.

## Features

- **Pangenome-based simulation**: Works with diploid samples from pangenome FASTA files
- **Realistic aDNA characteristics**: 
  - Configurable DNA fragment lengths (typical aDNA: ~45bp)
  - Single or double-stranded deamination patterns
  - C-to-T transitions at fragment ends
- **Contamination modeling**: Simulate exogenous human DNA contamination from other samples
- **Flexible sequencing**: Single-end (SE) or paired-end (PE) library support

## Installation

### Clone

```bash
git clone --recursive https://github.com/davidebolo1993/ancestralsim.git
cd ancestralsim
```

### Using conda/mamba (recommended)

```bash
ENV_PATH="/path/to/environment/installation/directory"
mamba env create -f environment.yml -p $ENV_PATH
#fix incorrext lib
cd $ENV_PATH/lib
ln -s libgsl.so.27 libgsl.so.25
cd -
```

### Using singularity

```bash
singularity pull --dir gargamel/. docker://quay.io/biocontainers/gargammel:1.1.4--hb66fcc3_0
```

When using this `singularity` container, you need to also have `bwa`, `fastp` and `samtools` executables available in `$PATH`.

## Quick Start

```bash
#Help
./ancestralsim.sh -h

# Basic usage
./ancestralsim.sh
  -r reference/GRCh38.fa
  -g ./gargammel
  -b region_mapping.tsv

# Increasing 0.5X coverage, increase modern-human contamination to 15%, simulate a single-end short-read library
./ancestralsim.sh
  -r reference/GRCh38.fa
  -g ./gargammel
  -b region_mapping.tsv
  -c 1.0
  --cont-ratio 0.15
  -o simulation_output
```

## Usage

```
Usage: ancestralsim.sh -r <reference.fa> -g <gargammel_dir> -b <region_mapping.tsv> [options]

Required:
-r FILE Reference genome FASTA
-g DIR Path to gargammel installation directory
-b FILE TSV mapping regions to alleles (chr<TAB>region_name<TAB>fasta_path)

Optional:
-c FLOAT Target coverage (default 0.5)
-l INT Mean fragment length (default 70)
-R INT Read length (default 100)
-L TYPE Library type: se|pe (default pe)
-d TYPE Deamination type: single|double (default single)
--deam-rate "VALS" Custom deam rates like "0.03,0.4,0.01,0.3"
--cont-ratio FLOAT Exogenous DNA ratio (default 0.1)
-o DIR Output directory (default output)
-t INT Threads (default 4)
-h Show this help
```

## Input Format

### Region Mapping File

Tab-separated file specifying chromosome, region name, and path to alleles FASTA:

```tsv
chr1  CYP2J2  /path/to/chr1_region1.fasta.gz
chr1  MUC1  /path/to/chr1_region2.fasta.gz
chr17 BRCA1 /path/to/chr17_region1.fasta.gz
chr22 CYP2D6  /path/to/chr22_region1.fasta.gz
```

### Alleles FASTA Format

The input FASTA file must contain diploid samples with headers in the format:

```
>SAMPLE_ID#HAPLOTYPE#CHROMOSOME
```

Example:
```
>HG00096#1#chr1
ACGTACGTACGT...
>HG00096#2#chr1
ACGTACGTACGT...
>HG00097#1#chr1
ACGTACGTACGT...
>HG00097#2#chr1
ACGTACGTACGT...
```

The pipeline automatically identifies diploid samples (those with both haplotype #1 and #2) and while it uses a diploid sample for simulation of aDNA reads, another random diploid sample from the same set is used to simulate modern-day contaminants. 

## Output

The pipeline generates:

```
output/
├── chr1/
│ ├── reference_chr1.fa # Extracted chromosome reference
│ ├── CYP2J2/
│ │ ├── bams/
│ │ │ ├── HG00096.sorted.bam
│ │ │ ├── HG00096.sorted.bam.bai
│ │ │ └── HG00096.flagstat.txt
│ │ ├── logs/
│ │ │ ├── HG00096_sequence_mapping.txt
│ │ │ ├── HG00096_gargammel.log
│ │ │ └── HG00096_fastp.html
│ │ └── region_summary.txt
│ └── MUC1/
│ └── (same structure)
├── chr17/
│ └── BRCA1/
│ └── (same structure)
├── merged_bams/ # Created after running merge script
│ ├── HG00096.merged.bam # Reads from all regions
│ └── HG00097.merged.bam
├── merge_commands.sh # Auto-generated merge script
└── global_simulation_summary.txt
```

### Key Output Files

- **`chr*/region/bams/*.sorted.bam`**: Per-region aligned aDNA reads
- **`merged_bams/*.merged.bam`**: Per-sample BAMs merged across all regions
- **`logs/*_sequence_mapping.txt`**: Maps simulated sequences to original haplotype names
- **`logs/*_fastp.html`**: Adapter trimming quality reports
- **`merge_commands.sh`**: Script to merge per-region BAMs into per-sample BAMs

## Merging BAMs Across Regions

After simulation completes, merge per-region BAMs into single BAM files per sample:

```
bash output/merge_commands.sh
```

This creates `output/merged_bams/SAMPLE.merged.bam` files containing reads from all simulated regions.


## Parameters Guide

### Coverage (`-c`)
Target sequencing coverage. For aDNA studies, typical values range from 0.1x to 2x.

### Fragment Length (`-l`)
Mean DNA fragment length. Ancient DNA is typically highly fragmented:
- **Modern DNA**: 150-500bp
- **Ancient DNA**: 30-80bp (mean ~45bp)

### Deamination (`-d`)
DNA damage patterns:
- **`single`**: C-to-T deamination on one strand (partially treated with UDG)
- **`double`**: C-to-T on both strands (non-UDG treated)

### Contamination (`--cont-ratio`)
Proportion of exogenous human DNA contamination (0.0-1.0):
- **0.0**: No contamination
- **0.1**: 10% contamination (typical for well-preserved samples)
- **0.3**: 30% contamination (moderate)
- **0.5+**: High contamination (challenging samples)

As mentioned above, the pipeline randomly selects a different diploid sample as the contaminant source.

## Citation

If you use ancestralsim in your research, please cite the gargammel paper:

- **gargammel**: Renaud G, Hanghøj K, Korneliussen TS, et al. (2017) gargammel: a sequence simulator for ancient DNA. *Bioinformatics* 33(4):577-579. [doi:10.1093/bioinformatics/btw670](https://doi.org/10.1093/bioinformatics/btw670)

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

For issues, questions, or suggestions, please open an issue on [GitHub](https://github.com/davidebolo1993/ancestralsim/issues).

## Author

Davide Bolognini ([@davidebolo1993](https://github.com/davidebolo1993))
