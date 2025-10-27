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

When using this `singularity` container, you need to also have `bwa` and `samtools` executables available in `$PATH`.

## Quick Start

```bash
#Help
./ancestralsim.sh -h

# Basic usage
./ancestralsim.sh \
  -f test/allele/chr22_39095117_39122647.fasta.gz \
  -r test/chr/chr22.fa.gz \
  -g ./gargammel \

# Increasing 0.5X coverage, increase modern-human contamination to 15%, simulate a single-end short-read library
./ancestralsim.sh \
  -f test/allele/chr22_39095117_39122647.fasta.gz \
  -r test/chr/chr22.fa.gz \
  -g ./gargammel \
  -c 1.0 \
  --cont-ratio 0.15 \
  -L se \
  -o test_output
```

## Usage

```
Usage: ancestralsim.sh -f <alleles.fasta> -r <reference.fa> -g <gargammel_dir> [options]

Required:
  -f FILE     Input FASTA file with pangenome alleles (can be gzipped)
  -r FILE     Reference chromosome FASTA for alignment
  -g DIR      Path to gargammel installation directory

Optional:
  -c FLOAT    Target coverage (default 0.5)
  -l INT      Mean fragment length (default 45)
  -R INT      Read length (default 75)
  -L TYPE     Library type: se|pe (default pe)
  -d TYPE     Deamination type: single|double (default single)
  --deam-rate "VALS"  Custom deam rates like "0.03,0.4,0.01,0.3"
  --cont-ratio FLOAT  Exogenous DNA ratio (default 0.1)
  -o DIR      Output directory (default output)
  -t INT      Threads (default 4)
  -h          Show this help
```

## Input Format

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
├── bams/
│   ├── SAMPLE1.sorted.bam
│   ├── SAMPLE1.sorted.bam.bai
│   ├── SAMPLE1.flagstat.txt
│   └── ...
├── logs/
│   ├── SAMPLE1_sequence_mapping.txt
│   ├── SAMPLE1_gargammel.log
│   ├── SAMPLE1_bwa_aln_1.log
│   └── ...
├── temp/
│   └── (intermediate files)
└── simulation_summary.txt
```

### Key Output Files

- **`bams/*.sorted.bam`**: Aligned aDNA reads in BAM format
- **`logs/*_sequence_mapping.txt`**: Maps simulated sequences back to original haplotype names
- **`logs/*_gargammel.log`**: Gargammel simulation logs
- **`simulation_summary.txt`**: Overall simulation statistics

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
