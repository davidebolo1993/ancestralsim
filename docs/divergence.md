# Extra Divergence

`--diverge` adds extra substitutions to selected endogenous haplotypes before gargammel simulates damaged reads. The selected PanSN haplotypes remain the background sequence; msprime is used to generate additional substitutions.

## Shared Options

```bash
--diverge
--divergence-rate FLOAT
--divergence-generations FLOAT
--divergence-years FLOAT
--generation-time FLOAT
--divergence-mode terminal|split
--divergence-model jc69
--seed INT
```

`--divergence-rate` is a mutation rate per bp per generation. A reasonable human nuclear DNA range is about `1.0e-8` to `1.4e-8`, with `1.25e-8` as the default.

You can specify time either directly in generations:

```bash
--divergence-generations 207
```

or in years:

```bash
--divergence-years 6000 --generation-time 29
```

These are mutually exclusive. The years form is converted internally:

```text
generations = divergence_years / generation_time
```

## Terminal Mode

`--divergence-mode terminal` is the default. It keeps the selected pangenome haplotypes and uses msprime to place mutations on one extra terminal branch per haplotype.

Expected substitutions are approximately:

```text
haplotype_length_bp * divergence_rate * sample_generations
```

For example, at `1.25e-8` and 29 years/generation, a 6000-year-old sample has about one expected added substitution per 386 kb. For short loci, biologically realistic divergence may often add zero mutations by chance.

## Split-Demography Mode

`--divergence-mode split` uses msprime's demography API to simulate one present-day haploid lineage and one ancient haploid lineage under a simple population split.

The sequence length passed to msprime is the length of each selected haplotype. For example, if `SAMPLE#1#chr1` is 120 kb long, msprime simulates a 120 kb ancestry for that haplotype.

The model is:

```text
ancestral population
        |
        | split time (--split-years or --split-generations)
       / \
modern   ancient
  |        |
today    ancient sample age (--divergence-years or --divergence-generations)
```

The ancestry is stochastic, but it is constrained by:

- ancient sample age,
- split time,
- effective population sizes,
- recombination rate,
- mutation rate.

msprime does not use the input FASTA as a reference genome. Instead, AncestralSim uses msprime to generate a difference pattern and overlays that pattern onto the selected real PanSN haplotype:

```text
1. Simulate one modern haploid lineage and one ancient haploid lineage.
2. Add mutations to that simulated ancestry with --divergence-rate.
3. Compare the simulated ancient and modern lineages.
4. Record positions where they differ.
5. Overlay those differences onto the selected real PanSN haplotype.
```

If msprime produces:

```text
position 120: simulated modern A, simulated ancient G
```

and the selected real haplotype has:

```text
position 120: C
```

AncestralSim changes the selected real haplotype at that position to:

```text
position 120: G
```

If the simulated ancient allele is already equal to the real base, AncestralSim chooses a different A/C/G/T base so that the overlay still creates an added difference.

## Split Population Sizes

By default, split mode uses `--effective-pop-size 10000` for all three populations:

```bash
--modern-ne 10000
--ancient-ne 10000
--ancestral-ne 10000
```

This is a conservative, common order-of-magnitude human coalescent effective population size. It is not a census population size.

For a simple sensitivity model with different branch sizes:

```bash
--modern-ne 10000 \
--ancient-ne 5000 \
--ancestral-ne 15000
```

This represents a standard modern branch, a smaller ancient branch with stronger drift, and a moderately larger ancestral population.

## Reports

Summary report:

```text
output/chr/region/logs/SAMPLE_divergence.tsv
```

Important columns:

- `mode`: `terminal` or `split`.
- `source`: `terminal_branch` or `split_demography`.
- `expected_mutations`: direct expectation in terminal mode; `NA` in split mode.
- `candidate_mutations`: simulated candidate differences.
- `applied_mutations`: substitutions overlaid onto the real haplotype.
- `modern_ne`, `ancient_ne`, `ancestral_ne`: split-mode population sizes.

Per-mutation report:

```text
output/chr/region/logs/SAMPLE_divergence_mutations.tsv
```

Important columns:

- `source`: `terminal_branch` or `split_demography`.
- `haplotype`: original PanSN haplotype name.
- `position_0based`, `position_1based`: modified position.
- `original_base`, `new_base`: actual modification applied to the selected haplotype.
- `simulated_modern`, `simulated_ancient`: alleles from the msprime simulation.

The diverged haplotypes with original PanSN headers are written to:

```text
output/chr/region/temp/SAMPLE/diverged_haplotypes.pansn.fa
```

