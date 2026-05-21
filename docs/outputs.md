# Outputs

The output directory is organized by chromosome, region, and sample.

## Layout

```text
output/
├── chr1/
│   ├── reference_chr1.fa
│   ├── CYP2J2/
│   │   ├── bams/
│   │   │   ├── HG00096.sorted.bam
│   │   │   ├── HG00096.sorted.bam.bai
│   │   │   └── HG00096.flagstat.txt
│   │   ├── logs/
│   │   │   ├── HG00096_sequence_mapping.txt
│   │   │   ├── HG00096_gargammel.log
│   │   │   ├── HG00096_fastp.html
│   │   │   ├── HG00096_divergence.tsv
│   │   │   └── HG00096_divergence_mutations.tsv
│   │   ├── temp/
│   │   └── region_summary.txt
│   └── neutral_control/
│       ├── neutral_control.pansn.fa
│       ├── bams/
│       ├── logs/
│       └── region_summary.txt
├── merged_bams/
│   ├── HG00096.merged.bam
│   └── HG00096.merged.bam.bai
├── merge_commands.sh
└── global_simulation_summary.txt
```

The divergence files are present only when `--diverge` is enabled. The control-region directory is present only when `--control-region` is provided.

## Key Files

Per-region BAMs:

```text
output/chr/region/bams/SAMPLE.sorted.bam
output/chr/region/bams/SAMPLE.sorted.bam.bai
output/chr/region/bams/SAMPLE.flagstat.txt
```

Per-sample logs:

```text
output/chr/region/logs/SAMPLE_sequence_mapping.txt
output/chr/region/logs/SAMPLE_gargammel.log
output/chr/region/logs/SAMPLE_fastp.html
output/chr/region/logs/SAMPLE_fastp.json
```

Divergence outputs:

```text
output/chr/region/logs/SAMPLE_divergence.tsv
output/chr/region/logs/SAMPLE_divergence_mutations.tsv
output/chr/region/temp/SAMPLE/diverged_haplotypes.pansn.fa
```

Control outputs:

```text
output/chr/control_name/control_name.pansn.fa
output/chr/control_name/bams/SAMPLE.sorted.bam
```

Run summary:

```text
output/global_simulation_summary.txt
```

## Merging BAMs

After simulation, merge per-region BAMs into one BAM per sample:

```bash
bash output/merge_commands.sh
```

Merged outputs:

```text
output/merged_bams/SAMPLE.merged.bam
output/merged_bams/SAMPLE.merged.bam.bai
```

The merge script includes all locus BAMs and, if requested, the control-region BAM.

