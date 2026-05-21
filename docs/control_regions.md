# Control Regions

Control regions simulate reads from a haploid reference interval and merge them with the locus BAMs for each sample.

## Usage

```bash
./ancestralsim.sh \
  -r reference/GRCh38.fa \
  -g ./gargammel \
  -b region_mapping.tsv \
  --control-region chr1:1000000-1050000 \
  --control-name neutral_control \
  -o simulation_output
```

`--control-region` must be formatted as:

```text
chr:start-end
```

`--control-name` controls the output directory name and defaults to `control`.

## What It Does

For every sample discovered in the locus FASTAs, the pipeline:

1. extracts the requested reference interval,
2. writes it as haploid endogenous gargammel input,
3. runs gargammel with the same coverage, fragment length, read length, library type, deamination, and contamination settings as the main loci,
4. aligns reads back to the extracted chromosome reference,
5. writes a per-sample control BAM,
6. adds the control BAM to `merge_commands.sh`.

Coverage is computed from the haploid reference interval in `endo/`. For example, `-c 0.5` means approximately 0.5x expected endogenous coverage over the control interval.

## Outputs

Control-region output is placed under:

```text
output/chr/control_name/
```

Important files:

```text
output/chr/control_name/control_name.pansn.fa
output/chr/control_name/bams/SAMPLE.sorted.bam
output/chr/control_name/bams/SAMPLE.sorted.bam.bai
output/chr/control_name/bams/SAMPLE.flagstat.txt
output/chr/control_name/logs/SAMPLE_sequence_mapping.txt
output/chr/control_name/region_summary.txt
```

The control FASTA header is written in PanSN-style form:

```text
>control_name#1#chr_start_end
```

