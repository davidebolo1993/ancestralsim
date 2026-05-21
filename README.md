# ancestralsim

Simulate ancient DNA (aDNA) reads from pangenome haplotypes.

## Overview

`ancestralsim` simulates aDNA sequencing reads from locus-specific diploid pangenome haplotypes. It uses [gargammel](https://github.com/grenaud/gargammel) for fragmentation, deamination, contamination, and Illumina read simulation, and can optionally use [msprime](https://tskit.dev/msprime/docs/stable/intro.html) to add extra terminal-branch or split-demography divergence before read simulation. Reads are trimmed with `fastp` and aligned to the reference chromosome with BWA.

Optional modules add:

- haploid reference control-region simulation,
- msprime-modeled extra divergence before read simulation,
- per-sample BAM merging across loci and controls.

## Documentation

Detailed documentation lives in [docs/](docs/README.md):

- [Inputs](docs/inputs.md): mapping file and PanSN FASTA requirements.
- [Core Simulation](docs/simulation.md): what happens for each locus/sample.
- [Control Regions](docs/control_regions.md): `--control-region` behavior and outputs.
- [Divergence](docs/divergence.md): terminal and split-demography divergence modes.
- [Outputs](docs/outputs.md): output layout and reports.
- [Example Scenarios](docs/examples.md): example commands and revision conditions.

## Installation

Clone recursively:

```bash
git clone --recursive https://github.com/davidebolo1993/ancestralsim.git
cd ancestralsim
```

Create the conda environment:

```bash
ENV_PATH="/path/to/environment/installation/directory"
mamba env create -f environment.yml -p "$ENV_PATH"
```

If your gargammel build needs the older GSL SONAME:

```bash
cd "$ENV_PATH/lib"
ln -s libgsl.so.27 libgsl.so.25 2>/dev/null || true
cd -
```

## Quick Start

Show help:

```bash
./ancestralsim.sh -h
```

Basic run:

```bash
./ancestralsim.sh \
  -r reference/GRCh38.fa \
  -g ./gargammel \
  -b region_mapping.tsv \
  -o simulation_output
```

Run with a control region and split-demography divergence:

```bash
./ancestralsim.sh \
  -r reference/GRCh38.fa \
  -g ./gargammel \
  -b region_mapping.tsv \
  --control-region chr1:1000000-1050000 \
  --control-name neutral_control \
  --diverge \
  --divergence-mode split \
  --divergence-years 5000 \
  --split-years 20000 \
  --generation-time 29 \
  --modern-ne 10000 \
  --ancient-ne 10000 \
  --ancestral-ne 10000 \
  -o simulation_output
```

Merge per-region BAMs after simulation:

```bash
bash simulation_output/merge_commands.sh
```

## Usage

```text
Usage: ancestralsim.sh -r <reference.fa> -g <gargammel_dir> -b <region_mapping.tsv> [options]
```

Required:

- `-r FILE`: reference genome FASTA.
- `-g DIR`: gargammel directory, used for matrices and local assets.
- `-b FILE`: region mapping TSV.

Common options:

- `-c FLOAT`: target endogenous coverage, default `0.5`.
- `-l INT`: mean fragment length, default `70`.
- `-R INT`: read length, default `100`.
- `-L se|pe`: library type, default `pe`.
- `-d single|double`: deamination profile, default `single`.
- `--cont-ratio FLOAT`: modern-human contamination ratio, default `0.1`.
- `--control-region chr:start-end`: add a haploid reference control interval.
- `--diverge`: enable extra divergence before read simulation.
- `--divergence-mode terminal|split`: divergence model, default `terminal`.
- `--divergence-years FLOAT`: ancient sample age in years.
- `--split-years FLOAT`: population split time in years for split mode.
- `--modern-ne`, `--ancient-ne`, `--ancestral-ne`: split-mode effective population sizes.
- `-o DIR`: output directory, default `output`.
- `-t INT`: threads, default `4`.

Run `./ancestralsim.sh -h` for the full option list.

## Minimal Input Example

Region mapping TSV:

```tsv
chr1    CYP2J2  /path/to/chr1_CYP2J2.fasta.gz
chr17   BRCA1   /path/to/chr17_BRCA1.fasta.gz
```

Haplotype FASTA headers:

```text
>HG00096#1#chr1
ACGT...
>HG00096#2#chr1
ACGT...
```

See [Inputs](docs/inputs.md) for full details.

## Citation

If you use ancestralsim in your research, please cite the relevant simulator papers:

- **gargammel**: Renaud G, Hanghøj K, Korneliussen TS, et al. (2017) gargammel: a sequence simulator for ancient DNA. *Bioinformatics* 33(4):577-579. [doi:10.1093/bioinformatics/btw670](https://doi.org/10.1093/bioinformatics/btw670)
- **msprime**: Kelleher J, Etheridge AM, McVean G. (2016) Efficient coalescent simulation and genealogical analysis for large sample sizes. *PLOS Computational Biology* 12(5):e1004842. [doi:10.1371/journal.pcbi.1004842](https://doi.org/10.1371/journal.pcbi.1004842)

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE).

## Author

Davide Bolognini ([@davidebolo1993](https://github.com/davidebolo1993))
