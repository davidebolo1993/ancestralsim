# Core Simulation

For each row in the region mapping file, AncestralSim simulates aDNA reads for every diploid sample in the locus FASTA.

## Per-Locus Flow

1. Extract the target chromosome from the reference FASTA.
2. Index the chromosome reference with BWA and samtools.
3. Index the locus haplotype FASTA if needed.
4. Identify diploid samples from PanSN-style FASTA headers.
5. For each sample:
   - extract haplotype `#1` and `#2` as endogenous input,
   - optionally apply extra divergence,
   - choose a different diploid sample as contaminant when possible,
   - create gargammel `endo/`, `cont/`, and `bact/` input directories,
   - run gargammel,
   - trim reads with fastp,
   - align reads to the extracted chromosome with BWA,
   - sort, index, and summarize the BAM with samtools.

## Endogenous Haplotype Handling

The selected endogenous haplotypes are first written with their original PanSN headers under:

```text
output/chr/region/temp/SAMPLE/source/
```

For gargammel, the pipeline writes copies with a shared internal header:

```text
>chr_endo
```

This is needed because gargammel expects the two endogenous chromosome FASTAs to contain matching sequence names.

## Contamination

When `--cont-ratio` is greater than zero, a different diploid sample from the same locus FASTA is randomly selected as the modern-human contaminant. If there is no other diploid sample, the target sample is reused.

The gargammel composition is:

```text
--comp 0,cont_ratio,1-cont_ratio
```

This means no bacterial contamination, the requested modern-human contamination fraction, and the remaining fraction as endogenous material.

## Main Parameters

- `-c`: endogenous target coverage.
- `-l`: mean fragment length.
- `-R`: read length.
- `-L`: library type, `pe` or `se`.
- `-d`: deamination profile, `single` or `double`.
- `--deam-rate`: custom Briggs damage parameters.
- `--cont-ratio`: modern-human contaminant fraction.

## Main Outputs

Per-region BAMs:

```text
output/chr/region/bams/SAMPLE.sorted.bam
output/chr/region/bams/SAMPLE.sorted.bam.bai
output/chr/region/bams/SAMPLE.flagstat.txt
```

Per-region logs:

```text
output/chr/region/logs/SAMPLE_sequence_mapping.txt
output/chr/region/logs/SAMPLE_gargammel.log
output/chr/region/logs/SAMPLE_fastp.html
```

See [Outputs](outputs.md) for the full layout.

