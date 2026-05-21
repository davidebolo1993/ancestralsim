# Inputs

AncestralSim requires a reference genome, a gargammel installation path, and a TSV that maps loci to haplotype FASTA files.

## Region Mapping File

The mapping file is tab-separated:

```tsv
chromosome  region_name  fasta_path
```

Example:

```tsv
chr1    CYP2J2  /path/to/chr1_CYP2J2.fasta.gz
chr1    MUC1    /path/to/chr1_MUC1.fasta.gz
chr17   BRCA1   /path/to/chr17_BRCA1.fasta.gz
chr22   CYP2D6  /path/to/chr22_CYP2D6.fasta.gz
```

The first field must match a sequence name in the reference FASTA.

## Haplotype FASTA

Each locus FASTA must contain diploid samples in PanSN-style headers:

```text
>SAMPLE_ID#HAPLOTYPE#CONTIG
```

Example:

```text
>HG00096#1#chr1
ACGTACGTACGT...
>HG00096#2#chr1
ACGTACGTACGT...
>HG00097#1#chr1
ACGTACGTACGT...
>HG00097#2#chr1
ACGTACGTACGT...
```

The pipeline identifies diploid samples as samples with exactly two records, `#1#` and `#2#`, in the FASTA index. When `--diverge` is used, these headers are validated and preserved in the diverged PanSN FASTA output.

## Required Tools

The recommended installation path is the conda environment in `environment.yml`. The runtime commands used by the wrapper include:

- `gargammel`
- `samtools`
- `bwa`
- `fastp`
- `python3`
- `bc`

