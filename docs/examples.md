# Example Scenarios

## Basic Run

```bash
./ancestralsim.sh \
  -r reference/GRCh38.fa \
  -g ./gargammel \
  -b region_mapping.tsv \
  -o simulation_output
```

## Higher Coverage and Contamination

```bash
./ancestralsim.sh \
  -r reference/GRCh38.fa \
  -g ./gargammel \
  -b region_mapping.tsv \
  -c 1.0 \
  --cont-ratio 0.15 \
  -o simulation_output
```

## Control Region

```bash
./ancestralsim.sh \
  -r reference/GRCh38.fa \
  -g ./gargammel \
  -b region_mapping.tsv \
  --control-region chr1:1000000-1050000 \
  --control-name neutral_control \
  -o simulation_output
```

## Terminal Divergence

```bash
./ancestralsim.sh \
  -r reference/GRCh38.fa \
  -g ./gargammel \
  -b region_mapping.tsv \
  --diverge \
  --divergence-mode terminal \
  --divergence-years 6000 \
  --generation-time 29 \
  --seed 13 \
  -o terminal_6k
```

## Split-Demography Revision Conditions

Conservative baseline with equal effective population sizes:

```bash
# anc5k_split20k
./ancestralsim.sh \
  -r reference/GRCh38.fa \
  -g ./gargammel \
  -b region_mapping.tsv \
  --diverge \
  --divergence-mode split \
  --divergence-years 5000 \
  --split-years 20000 \
  --generation-time 29 \
  --modern-ne 10000 \
  --ancient-ne 10000 \
  --ancestral-ne 10000 \
  -o anc5k_split20k

# anc20k_split100k
./ancestralsim.sh \
  -r reference/GRCh38.fa \
  -g ./gargammel \
  -b region_mapping.tsv \
  --diverge \
  --divergence-mode split \
  --divergence-years 20000 \
  --split-years 100000 \
  --generation-time 29 \
  --modern-ne 10000 \
  --ancient-ne 10000 \
  --ancestral-ne 10000 \
  -o anc20k_split100k
```

Exploratory sensitivity setting with different effective population sizes:

```bash
--modern-ne 10000 \
--ancient-ne 5000 \
--ancestral-ne 15000
```

Use this when you explicitly want a smaller ancient branch and a larger ancestral branch.
