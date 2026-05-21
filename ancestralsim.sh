#!/bin/bash
# ancestralsim: Ancient DNA simulation from pangenome haplotypes

# Defaults
COVERAGE=0.5
FRAGMENT_LENGTH=70
READ_LENGTH=100
DEAMINATION="single"
DEAM_RATE=""
LIBRARY_TYPE="pe"
REFERENCE=""
OUTPUT_DIR="output"
THREADS=4
GARGAMMEL_DIR=""
CONT_RATIO=0.1
REGION_MAP=""
CONTROL_REGION=""
CONTROL_NAME="control"
DIVERGE=false
DIVERGENCE_RATE="1.25e-8"
DIVERGENCE_GENERATIONS=1000
DIVERGENCE_GENERATIONS_SET=false
DIVERGENCE_YEARS=""
GENERATION_TIME=29
DIVERGENCE_MODEL="jc69"
DIVERGENCE_MODE="terminal"
SPLIT_GENERATIONS=""
SPLIT_YEARS=""
EFFECTIVE_POP_SIZE=10000
MODERN_NE=""
ANCIENT_NE=""
ANCESTRAL_NE=""
RECOMBINATION_RATE="1e-8"
SEED=""

usage() {
local exit_code="${1:-1}"
cat <<EOF
Usage: $0 -r <reference.fa> -g <gargammel_dir> -b <region_mapping.tsv> [options]

Required:
  -r FILE     Reference genome FASTA
  -g DIR      Path to gargammel installation directory
  -b FILE     TSV mapping regions to alleles (chr<TAB>region_name<TAB>fasta_path)

Optional:
  -c FLOAT    Target coverage (default 0.5)
  -l INT      Mean fragment length (default 70)
  -R INT      Read length (default 100)
  -L TYPE     Library type: se|pe (default pe)
  -d TYPE     Deamination type: single|double (default single)
  --deam-rate "VALS"  Custom deam rates like "0.03,0.4,0.01,0.3"
  --cont-ratio FLOAT  Exogenous DNA ratio (default 0.1)
  --control-region REGION  Haploid reference control region chr:start-end
  --control-name NAME      Name for control-region outputs (default control)
  --diverge               Add msprime-modeled divergence to endogenous haplotypes
  --divergence-rate FLOAT Mutation rate per bp per generation (default 1.25e-8)
  --divergence-generations INT  Terminal branch length in generations (default 1000)
  --divergence-years FLOAT      Sample age/divergence time in years before present
  --generation-time FLOAT       Years per generation for --divergence-years (default 29)
  --divergence-mode MODE        terminal|split (default terminal)
  --split-generations FLOAT     Population split time in generations (split mode)
  --split-years FLOAT           Population split time in years before present (split mode)
  --effective-pop-size FLOAT    Set all split-mode Ne values (default 10000)
  --modern-ne FLOAT             Modern population Ne for split mode (default: --effective-pop-size)
  --ancient-ne FLOAT            Ancient population Ne for split mode (default: --effective-pop-size)
  --ancestral-ne FLOAT          Ancestral population Ne for split mode (default: --effective-pop-size)
  --recombination-rate FLOAT    Recombination rate for split mode (default 1e-8)
  --divergence-model MODEL     Mutation model: jc69 (default jc69)
  --seed INT             Seed for reproducible divergence
  -o DIR      Output directory (default output)
  -t INT      Threads (default 4)
  -h          Show this help

EOF
exit "$exit_code"
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -r) REFERENCE="$2"; shift 2 ;;
        -g) GARGAMMEL_DIR="$2"; shift 2 ;;
        -b) REGION_MAP="$2"; shift 2 ;;
        -c) COVERAGE="$2"; shift 2 ;;
        -l) FRAGMENT_LENGTH="$2"; shift 2 ;;
        -R) READ_LENGTH="$2"; shift 2 ;;
        -L) LIBRARY_TYPE="$2"; shift 2 ;;
        -d) DEAMINATION="$2"; shift 2 ;;
        --deam-rate) DEAM_RATE="$2"; shift 2 ;;
        --cont-ratio) CONT_RATIO="$2"; shift 2 ;;
        --control-region) CONTROL_REGION="$2"; shift 2 ;;
        --control-name) CONTROL_NAME="$2"; shift 2 ;;
        --diverge) DIVERGE=true; shift ;;
        --divergence-rate) DIVERGENCE_RATE="$2"; shift 2 ;;
        --divergence-generations) DIVERGENCE_GENERATIONS="$2"; DIVERGENCE_GENERATIONS_SET=true; shift 2 ;;
        --divergence-years) DIVERGENCE_YEARS="$2"; shift 2 ;;
        --generation-time) GENERATION_TIME="$2"; shift 2 ;;
        --divergence-mode) DIVERGENCE_MODE="$2"; shift 2 ;;
        --split-generations) SPLIT_GENERATIONS="$2"; shift 2 ;;
        --split-years) SPLIT_YEARS="$2"; shift 2 ;;
        --effective-pop-size) EFFECTIVE_POP_SIZE="$2"; shift 2 ;;
        --modern-ne) MODERN_NE="$2"; shift 2 ;;
        --ancient-ne) ANCIENT_NE="$2"; shift 2 ;;
        --ancestral-ne) ANCESTRAL_NE="$2"; shift 2 ;;
        --recombination-rate) RECOMBINATION_RATE="$2"; shift 2 ;;
        --divergence-model) DIVERGENCE_MODEL="$2"; shift 2 ;;
        --seed) SEED="$2"; shift 2 ;;
        -o) OUTPUT_DIR="$2"; shift 2 ;;
        -t) THREADS="$2"; shift 2 ;;
        -h) usage 0 ;;
        *) echo "Unknown option $1"; usage ;;
    esac
done

if [[ -z "$REFERENCE" || -z "$GARGAMMEL_DIR" || -z "$REGION_MAP" ]]; then
    echo "Error: -r, -g, and -b are required."; usage
fi
if [[ ! -d "$GARGAMMEL_DIR" ]]; then
    echo "Error: Gargammel directory not found: $GARGAMMEL_DIR"; exit 1
fi
if [[ ! -f "$REGION_MAP" ]]; then
    echo "Error: Region mapping file not found: $REGION_MAP"; exit 1
fi

REFERENCE=$(realpath "$REFERENCE")
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR=$(cd "$OUTPUT_DIR" && pwd -P)
GARGAMMEL_DIR=$(realpath "$GARGAMMEL_DIR")
REGION_MAP=$(realpath "$REGION_MAP")
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

if command -v gargammel &>/dev/null; then
    GARGAMMEL_CMD="gargammel"
    echo "Using gargammel from PATH"
elif [[ -f "$GARGAMMEL_DIR/gargammel.pl" ]]; then
    GARGAMMEL_CMD="perl \"$GARGAMMEL_DIR/gargammel.pl\""
    echo "Using gargammel.pl from $GARGAMMEL_DIR"
else
    echo "Error: gargammel not found in PATH and gargammel.pl not found in $GARGAMMEL_DIR"; exit 1
fi

MATRIX_DIR="$GARGAMMEL_DIR/src/matrices"
MATRIX_SINGLE="$MATRIX_DIR/single-"
MATRIX_DOUBLE="$MATRIX_DIR/double-"

if [[ ! -d "$MATRIX_DIR" ]]; then
    echo "Error: Matrices directory missing: $MATRIX_DIR"; exit 1
fi

if [[ "$LIBRARY_TYPE" != "se" ]] && [[ "$LIBRARY_TYPE" != "pe" ]]; then
    echo "Error: Library type must be 'se' or 'pe'"; exit 1
fi

if [[ "$DEAMINATION" != "single" ]] && [[ "$DEAMINATION" != "double" ]]; then
    echo "Error: Deamination type must be 'single' or 'double'"; exit 1
fi

if (( $(echo "$CONT_RATIO < 0" | bc -l) )) || (( $(echo "$CONT_RATIO > 1" | bc -l) )); then
    echo "Error: --cont-ratio must be between 0 and 1"; exit 1
fi

ENDO_RATIO=$(python3 -c "print(round(1.0 - $CONT_RATIO, 6))")
MODERN_NE="${MODERN_NE:-$EFFECTIVE_POP_SIZE}"
ANCIENT_NE="${ANCIENT_NE:-$EFFECTIVE_POP_SIZE}"
ANCESTRAL_NE="${ANCESTRAL_NE:-$EFFECTIVE_POP_SIZE}"

if [[ "$DIVERGENCE_MODEL" != "jc69" ]]; then
    echo "Error: --divergence-model currently supports only 'jc69'"; exit 1
fi

if [[ "$DIVERGENCE_MODE" != "terminal" ]] && [[ "$DIVERGENCE_MODE" != "split" ]]; then
    echo "Error: --divergence-mode must be 'terminal' or 'split'"; exit 1
fi

if [[ -n "$DIVERGENCE_YEARS" ]]; then
    if [[ "$DIVERGENCE_GENERATIONS_SET" == true ]]; then
        echo "Error: use either --divergence-generations or --divergence-years, not both"; exit 1
    fi
    DIVERGENCE_GENERATIONS=$(python3 - "$DIVERGENCE_YEARS" "$GENERATION_TIME" <<'PY'
import sys
years = float(sys.argv[1])
generation_time = float(sys.argv[2])
if years < 0:
    raise SystemExit("--divergence-years must be >= 0")
if generation_time <= 0:
    raise SystemExit("--generation-time must be > 0")
print(years / generation_time)
PY
)
fi

if [[ -n "$SPLIT_YEARS" ]]; then
    if [[ -n "$SPLIT_GENERATIONS" ]]; then
        echo "Error: use either --split-generations or --split-years, not both"; exit 1
    fi
    SPLIT_GENERATIONS=$(python3 - "$SPLIT_YEARS" "$GENERATION_TIME" <<'PY'
import sys
years = float(sys.argv[1])
generation_time = float(sys.argv[2])
if years < 0:
    raise SystemExit("--split-years must be >= 0")
if generation_time <= 0:
    raise SystemExit("--generation-time must be > 0")
print(years / generation_time)
PY
)
fi

if [[ "$DIVERGENCE_MODE" == "split" ]]; then
    if [[ -z "$SPLIT_GENERATIONS" ]]; then
        echo "Error: --divergence-mode split requires --split-generations or --split-years"; exit 1
    fi
    python3 - "$DIVERGENCE_GENERATIONS" "$SPLIT_GENERATIONS" "$MODERN_NE" "$ANCIENT_NE" "$ANCESTRAL_NE" "$RECOMBINATION_RATE" <<'PY'
import sys
sample_generations = float(sys.argv[1])
split_generations = float(sys.argv[2])
modern_ne = float(sys.argv[3])
ancient_ne = float(sys.argv[4])
ancestral_ne = float(sys.argv[5])
recombination_rate = float(sys.argv[6])
if split_generations <= sample_generations:
    raise SystemExit("--split-generations/--split-years must be older than the ancient sample age")
if modern_ne <= 0:
    raise SystemExit("--modern-ne must be > 0")
if ancient_ne <= 0:
    raise SystemExit("--ancient-ne must be > 0")
if ancestral_ne <= 0:
    raise SystemExit("--ancestral-ne must be > 0")
if recombination_rate < 0:
    raise SystemExit("--recombination-rate must be >= 0")
PY
    if [[ $? -ne 0 ]]; then
        exit 1
    fi
fi

if [[ -n "$CONTROL_REGION" ]]; then
    if [[ "$CONTROL_REGION" =~ ^([^:]+):([0-9]+)-([0-9]+)$ ]]; then
        CONTROL_CHR="${BASH_REMATCH[1]}"
        CONTROL_START="${BASH_REMATCH[2]}"
        CONTROL_END="${BASH_REMATCH[3]}"
        CONTROL_CONTIG="${CONTROL_CHR}_${CONTROL_START}_${CONTROL_END}"
        if (( CONTROL_START < 1 || CONTROL_END < CONTROL_START )); then
            echo "Error: invalid --control-region coordinates: $CONTROL_REGION"; exit 1
        fi
    else
        echo "Error: --control-region must be formatted as chr:start-end"; exit 1
    fi
fi

for cmd in samtools bwa fastp python3 bc; do
    if ! command -v $cmd &> /dev/null; then
        echo "Error: $cmd not found. Please install it."; exit 1
    fi
done

if [[ "$DIVERGE" == true ]]; then
    if [[ ! -f "$SCRIPT_DIR/scripts/diverge_haplotypes_msprime.py" ]]; then
        echo "Error: divergence helper missing: $SCRIPT_DIR/scripts/diverge_haplotypes_msprime.py"; exit 1
    fi
    python3 -c "import msprime, tskit" 2>/dev/null || {
        echo "Error: --diverge requires msprime and tskit in the active environment"; exit 1
    }
fi

derive_seed() {
    local label="$1"
    if [[ -z "$SEED" ]]; then
        echo ""
    else
        local hashed
        hashed=$(printf "%s" "${SEED}:${label}" | cksum | awk '{print $1}')
        echo $(( (hashed % 2147483646) + 1 ))
    fi
}

mark_extracted_chr() {
    printf '%s\n' "$1" >> "$CHR_EXTRACTED_LIST"
}

have_extracted_chr() {
    grep -Fxq "$1" "$CHR_EXTRACTED_LIST"
}

mkdir -p "$OUTPUT_DIR"

if [[ ! -f "${REFERENCE}.fai" ]]; then
    echo "Indexing reference with samtools..."
    samtools faidx "$REFERENCE"
fi

echo "=== AncestralSim: Ancient DNA Simulation Pipeline ==="
echo "Reference: $REFERENCE"
echo "Region mapping: $REGION_MAP"
echo "Coverage: ${COVERAGE}x, Fragment: ${FRAGMENT_LENGTH}bp, Read: ${READ_LENGTH}bp"
echo "Library: $LIBRARY_TYPE, Deamination: $DEAMINATION, Contamination: ${CONT_RATIO}"
if [[ -n "$CONTROL_REGION" ]]; then
    echo "Control region: $CONTROL_REGION ($CONTROL_NAME)"
fi
if [[ "$DIVERGE" == true ]]; then
    echo "Divergence: enabled, mode=${DIVERGENCE_MODE}, model=${DIVERGENCE_MODEL}, rate=${DIVERGENCE_RATE}, generations=${DIVERGENCE_GENERATIONS}"
    if [[ -n "$DIVERGENCE_YEARS" ]]; then
        echo "Divergence years: ${DIVERGENCE_YEARS}, generation time: ${GENERATION_TIME}"
    fi
    if [[ "$DIVERGENCE_MODE" == "split" ]]; then
        echo "Split: generations=${SPLIT_GENERATIONS}, modern Ne=${MODERN_NE}, ancient Ne=${ANCIENT_NE}, ancestral Ne=${ANCESTRAL_NE}, recombination=${RECOMBINATION_RATE}"
    fi
fi
echo ""

CHR_EXTRACTED_LIST="$OUTPUT_DIR/.chromosomes_extracted.list"
ALL_SAMPLES_FILE="$OUTPUT_DIR/.all_samples.list"
UNIQUE_SAMPLES_FILE="$OUTPUT_DIR/.unique_samples.list"
: > "$CHR_EXTRACTED_LIST"
: > "$ALL_SAMPLES_FILE"
: > "$UNIQUE_SAMPLES_FILE"

N_REGIONS=$(grep -v "^#" "$REGION_MAP" | grep -v "^$" | wc -l | awk '{print $1}')
echo "Found $N_REGIONS regions to process"
echo ""

GLOBAL_SUMMARY="$OUTPUT_DIR/global_simulation_summary.txt"
MERGE_SCRIPT="$OUTPUT_DIR/merge_commands.sh"

cat > "$GLOBAL_SUMMARY" <<EOF
AncestralSim Simulation Report
===============================
Date: $(date)
Reference: $REFERENCE
Region mapping: $REGION_MAP

Parameters:
- Coverage: ${COVERAGE}x
- Fragment length: ${FRAGMENT_LENGTH}bp
- Read length: ${READ_LENGTH}bp
- Library: $LIBRARY_TYPE
- Deamination: $DEAMINATION
- Contamination: ${CONT_RATIO}
- Control region: ${CONTROL_REGION:-none}
- Divergence: $DIVERGE
- Divergence mode: $DIVERGENCE_MODE
- Divergence model: $DIVERGENCE_MODEL
- Divergence rate: $DIVERGENCE_RATE
- Divergence generations: $DIVERGENCE_GENERATIONS
- Divergence years: ${DIVERGENCE_YEARS:-none}
- Generation time: $GENERATION_TIME
- Split generations: ${SPLIT_GENERATIONS:-none}
- Split years: ${SPLIT_YEARS:-none}
- Effective population size: $EFFECTIVE_POP_SIZE
- Modern Ne: $MODERN_NE
- Ancient Ne: $ANCIENT_NE
- Ancestral Ne: $ANCESTRAL_NE
- Recombination rate: $RECOMBINATION_RATE
- Regions: $N_REGIONS

Per-Region Details:
====================

EOF

cat > "$MERGE_SCRIPT" <<'EOF'
#!/bin/bash
# Merge BAMs per sample across all regions

set -e

OUTPUT_DIR="$(dirname "$0")"
MERGED_DIR="$OUTPUT_DIR/merged_bams"
mkdir -p "$MERGED_DIR"

echo "=== Merging BAMs per sample ==="
echo "Output: $MERGED_DIR"
echo ""

EOF

REGION_NUM=0
while IFS=$'\t' read -r TARGET_CHR REGION_NAME ALLELES_FASTA; do
    [[ "$TARGET_CHR" =~ ^#.*$ ]] && continue
    [[ -z "$TARGET_CHR" || -z "$REGION_NAME" || -z "$ALLELES_FASTA" ]] && continue
    
    REGION_NUM=$((REGION_NUM + 1))
    
    echo "Processing region $REGION_NUM/$N_REGIONS: ${TARGET_CHR}/${REGION_NAME}"
    
    ALLELES_FASTA=$(realpath "$ALLELES_FASTA")
    
    if [[ ! -f "$ALLELES_FASTA" ]]; then
        echo "WARNING: File not found: $ALLELES_FASTA, skipping"
        echo ""
        continue
    fi
    
    REGION_OUTPUT="$OUTPUT_DIR/${TARGET_CHR}/${REGION_NAME}"
    mkdir -p "$REGION_OUTPUT"/{logs,temp,bams}
    
    CHROM_REFERENCE="$OUTPUT_DIR/${TARGET_CHR}/reference_${TARGET_CHR}.fa"
    
    if ! have_extracted_chr "$TARGET_CHR"; then
        echo "Extracting chromosome ${TARGET_CHR}..."
        mkdir -p "$OUTPUT_DIR/${TARGET_CHR}"
        samtools faidx "$REFERENCE" "$TARGET_CHR" > "$CHROM_REFERENCE"
        
        if [[ ! -s "$CHROM_REFERENCE" ]]; then
            echo "ERROR: Failed to extract $TARGET_CHR, skipping"
            echo ""
            continue
        fi
        
        echo "Indexing ${TARGET_CHR} with BWA..."
        bwa index "$CHROM_REFERENCE"
        samtools faidx "$CHROM_REFERENCE"
        
        mark_extracted_chr "$TARGET_CHR"
    fi
    
    ALIGN_REFERENCE="$CHROM_REFERENCE"
    
    if [[ ! -f "${ALLELES_FASTA}.fai" ]]; then
        samtools faidx "$ALLELES_FASTA"
    fi
    FAI_FILE="${ALLELES_FASTA}.fai"
    
    # Identify diploid samples (2 alleles per sample)
    awk -F'[#\t]' '{print $1}' "$FAI_FILE" | sort | uniq -c | \
    awk '$1==2 {print $2}' > "$REGION_OUTPUT/temp/diploid_samples.txt"
    NDIPLOID=$(wc -l < "$REGION_OUTPUT/temp/diploid_samples.txt" | awk '{print $1}')
    
    if [[ $NDIPLOID -eq 0 ]]; then
        echo "WARNING: No diploid samples found, skipping"
        echo ""
        continue
    fi
    
    echo "Found ${NDIPLOID} diploid samples"
    echo ""
    
    while IFS= read -r sample; do
        printf '%s\n' "$sample" >> "$ALL_SAMPLES_FILE"
    done < "$REGION_OUTPUT/temp/diploid_samples.txt"
    
    cat >> "$GLOBAL_SUMMARY" <<EOF
${TARGET_CHR}/${REGION_NAME}: $NDIPLOID samples
  Alleles: $ALLELES_FASTA
  Output: $REGION_OUTPUT/bams/

EOF
    
    while IFS= read -r sample; do
        echo "  Sample: ${sample}"
        SAMPLE_DIR="$REGION_OUTPUT/temp/${sample}"
        mkdir -p "$SAMPLE_DIR"/{endo,cont,bact,source}
        
        MAPPING_LOG="$REGION_OUTPUT/logs/${sample}_sequence_mapping.txt"
        echo "Sample: ${sample} | Chromosome: ${TARGET_CHR} | Region: ${REGION_NAME}" > "$MAPPING_LOG"
        echo "Date: $(date)" >> "$MAPPING_LOG"
        echo "" >> "$MAPPING_LOG"

        ENDO_HAP1_ORIG=$(grep "^${sample}#1#" "$FAI_FILE" | cut -f1)
        ENDO_HAP2_ORIG=$(grep "^${sample}#2#" "$FAI_FILE" | cut -f1)
        
        samtools faidx "$ALLELES_FASTA" "$ENDO_HAP1_ORIG" > "$SAMPLE_DIR/source/hap1.pansn.fa"
        samtools faidx "$ALLELES_FASTA" "$ENDO_HAP2_ORIG" > "$SAMPLE_DIR/source/hap2.pansn.fa"

        if [[ "$DIVERGE" == true ]]; then
            DIVERGENCE_SEED=$(derive_seed "${TARGET_CHR}:${REGION_NAME}:${sample}")
            DIVERGENCE_CMD=(python3 "$SCRIPT_DIR/scripts/diverge_haplotypes_msprime.py"
                --hap1 "$SAMPLE_DIR/source/hap1.pansn.fa"
                --hap2 "$SAMPLE_DIR/source/hap2.pansn.fa"
                --out-hap1 "$SAMPLE_DIR/endo/hap1.fa"
                --out-hap2 "$SAMPLE_DIR/endo/hap2.fa"
                --pansn-output "$SAMPLE_DIR/diverged_haplotypes.pansn.fa"
                --report "$REGION_OUTPUT/logs/${sample}_divergence.tsv"
                --mutation-report "$REGION_OUTPUT/logs/${sample}_divergence_mutations.tsv"
                --header chr_endo
                --mode "$DIVERGENCE_MODE"
                --rate "$DIVERGENCE_RATE"
                --generations "$DIVERGENCE_GENERATIONS"
                --effective-pop-size "$EFFECTIVE_POP_SIZE"
                --modern-ne "$MODERN_NE"
                --ancient-ne "$ANCIENT_NE"
                --ancestral-ne "$ANCESTRAL_NE"
                --recombination-rate "$RECOMBINATION_RATE"
                --model "$DIVERGENCE_MODEL"
                --require-pansn)
            if [[ "$DIVERGENCE_MODE" == "split" ]]; then
                DIVERGENCE_CMD+=(--split-generations "$SPLIT_GENERATIONS")
            fi
            if [[ -n "$DIVERGENCE_SEED" ]]; then
                DIVERGENCE_CMD+=(--seed "$DIVERGENCE_SEED")
            fi
            "${DIVERGENCE_CMD[@]}"
        else
            sed "s/^>.*/>chr_endo/" "$SAMPLE_DIR/source/hap1.pansn.fa" > "$SAMPLE_DIR/endo/hap1.fa"
            sed "s/^>.*/>chr_endo/" "$SAMPLE_DIR/source/hap2.pansn.fa" > "$SAMPLE_DIR/endo/hap2.fa"
        fi
        
        echo "[ENDOGENOUS]" >> "$MAPPING_LOG"
        echo "hap1.fa <- $ENDO_HAP1_ORIG" >> "$MAPPING_LOG"
        echo "hap2.fa <- $ENDO_HAP2_ORIG" >> "$MAPPING_LOG"
        if [[ "$DIVERGE" == true ]]; then
            echo "divergence <- msprime mode=${DIVERGENCE_MODE}, model=${DIVERGENCE_MODEL}, rate=${DIVERGENCE_RATE}, generations=${DIVERGENCE_GENERATIONS}" >> "$MAPPING_LOG"
            if [[ "$DIVERGENCE_MODE" == "split" ]]; then
                echo "split <- generations=${SPLIT_GENERATIONS}, modern_ne=${MODERN_NE}, ancient_ne=${ANCIENT_NE}, ancestral_ne=${ANCESTRAL_NE}, recombination=${RECOMBINATION_RATE}" >> "$MAPPING_LOG"
            fi
            echo "diverged PanSN FASTA: $SAMPLE_DIR/diverged_haplotypes.pansn.fa" >> "$MAPPING_LOG"
        fi
        echo "" >> "$MAPPING_LOG"

        if (( $(echo "$CONT_RATIO > 0" | bc -l) )); then
            OTHER_SAMPLES=$(grep -v "^${sample}$" "$REGION_OUTPUT/temp/diploid_samples.txt")
            
            if [[ -n "$OTHER_SAMPLES" ]]; then
                CONT_SAMPLE=$(echo "$OTHER_SAMPLES" | shuf -n 1)
            else
                CONT_SAMPLE="$sample"
            fi
            
            CONT_HAP1_ORIG=$(grep "^${CONT_SAMPLE}#1#" "$FAI_FILE" | cut -f1)
            CONT_HAP2_ORIG=$(grep "^${CONT_SAMPLE}#2#" "$FAI_FILE" | cut -f1)
            
            samtools faidx "$ALLELES_FASTA" "$CONT_HAP1_ORIG" | \
                sed "s/^>.*/>chr_cont/" > "$SAMPLE_DIR/cont/cont_hap1.fa"
            samtools faidx "$ALLELES_FASTA" "$CONT_HAP2_ORIG" | \
                sed "s/^>.*/>chr_cont/" > "$SAMPLE_DIR/cont/cont_hap2.fa"
            
            echo "[CONTAMINANT] - ${CONT_SAMPLE}" >> "$MAPPING_LOG"
            echo "cont_hap1.fa <- $CONT_HAP1_ORIG" >> "$MAPPING_LOG"
            echo "cont_hap2.fa <- $CONT_HAP2_ORIG" >> "$MAPPING_LOG"
            echo "" >> "$MAPPING_LOG"
        else
            samtools faidx "$ALLELES_FASTA" "$ENDO_HAP1_ORIG" | \
                sed "s/^>.*/>chr_cont/" > "$SAMPLE_DIR/cont/dummy.fa"
            echo "[CONTAMINANT] - None" >> "$MAPPING_LOG"
            echo "" >> "$MAPPING_LOG"
        fi

        samtools faidx "$ALLELES_FASTA" "$ENDO_HAP1_ORIG" | \
            sed "s/^>.*/>chr_bact/" > "$SAMPLE_DIR/bact/dummy.fa"
        echo "[BACTERIAL] - None" >> "$MAPPING_LOG"
        echo "" >> "$MAPPING_LOG"

        GARGAMMEL_RUN="$GARGAMMEL_CMD \
            -c $COVERAGE \
            --comp 0,$CONT_RATIO,$ENDO_RATIO \
            -l $FRAGMENT_LENGTH \
            -rl $READ_LENGTH \
            -o $SAMPLE_DIR/sim"
        
        if [[ -n "$DEAM_RATE" ]]; then
            GARGAMMEL_RUN="$GARGAMMEL_RUN -damage $DEAM_RATE"
        else
            if [[ "$DEAMINATION" == "single" ]]; then
                MATRIX="$MATRIX_SINGLE"
            else
                MATRIX="$MATRIX_DOUBLE"
            fi
            GARGAMMEL_RUN="$GARGAMMEL_RUN -matfile $MATRIX"
        fi
        
        if [[ "$LIBRARY_TYPE" == "se" ]]; then
            GARGAMMEL_RUN="$GARGAMMEL_RUN -se"
        fi
        
        GARGAMMEL_RUN="$GARGAMMEL_RUN $SAMPLE_DIR"
        eval $GARGAMMEL_RUN > "$REGION_OUTPUT/logs/${sample}_gargammel.log" 2>&1

        if [[ "$LIBRARY_TYPE" == "pe" ]]; then
            fastp \
                --in1 "$SAMPLE_DIR/sim_s1.fq.gz" \
                --in2 "$SAMPLE_DIR/sim_s2.fq.gz" \
                --out1 "$SAMPLE_DIR/sim_s1_trimmed.fq.gz" \
                --out2 "$SAMPLE_DIR/sim_s2_trimmed.fq.gz" \
                --detect_adapter_for_pe \
                --thread "$THREADS" \
                --length_required 25 \
                --json "$REGION_OUTPUT/logs/${sample}_fastp.json" \
                --html "$REGION_OUTPUT/logs/${sample}_fastp.html" \
                2> "$REGION_OUTPUT/logs/${sample}_fastp.log"
            
            READ1="$SAMPLE_DIR/sim_s1_trimmed.fq.gz"
            READ2="$SAMPLE_DIR/sim_s2_trimmed.fq.gz"
        else
            fastp \
                --in1 "$SAMPLE_DIR/sim_s.fq.gz" \
                --out1 "$SAMPLE_DIR/sim_s_trimmed.fq.gz" \
                --thread "$THREADS" \
                --length_required 25 \
                --json "$REGION_OUTPUT/logs/${sample}_fastp.json" \
                --html "$REGION_OUTPUT/logs/${sample}_fastp.html" \
                --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAG \
                2> "$REGION_OUTPUT/logs/${sample}_fastp.log"
            
            READ1="$SAMPLE_DIR/sim_s_trimmed.fq.gz"
        fi

        if [[ "$LIBRARY_TYPE" == "pe" ]]; then
            bwa aln -l 16500 -n 0.01 -o 2 -t "$THREADS" \
                "$ALIGN_REFERENCE" "$READ1" \
                > "$SAMPLE_DIR/sim_1.sai" \
                2> "$REGION_OUTPUT/logs/${sample}_bwa_aln_1.log"
            
            bwa aln -l 16500 -n 0.01 -o 2 -t "$THREADS" \
                "$ALIGN_REFERENCE" "$READ2" \
                > "$SAMPLE_DIR/sim_2.sai" \
                2> "$REGION_OUTPUT/logs/${sample}_bwa_aln_2.log"
            
            bwa sampe "$ALIGN_REFERENCE" \
                "$SAMPLE_DIR/sim_1.sai" "$SAMPLE_DIR/sim_2.sai" \
                "$READ1" "$READ2" \
                2> "$REGION_OUTPUT/logs/${sample}_bwa_sampe.log" | \
            samtools sort -@ "$THREADS" \
                -o "$REGION_OUTPUT/bams/${sample}.sorted.bam" -
        else
            bwa aln -l 16500 -n 0.01 -o 2 -t "$THREADS" \
                "$ALIGN_REFERENCE" "$READ1" \
                > "$SAMPLE_DIR/sim.sai" \
                2> "$REGION_OUTPUT/logs/${sample}_bwa_aln.log"
            
            bwa samse "$ALIGN_REFERENCE" \
                "$SAMPLE_DIR/sim.sai" "$READ1" \
                2> "$REGION_OUTPUT/logs/${sample}_bwa_samse.log" | \
            samtools sort -@ "$THREADS" \
                -o "$REGION_OUTPUT/bams/${sample}.sorted.bam" -
        fi

        samtools index "$REGION_OUTPUT/bams/${sample}.sorted.bam"
        samtools flagstat "$REGION_OUTPUT/bams/${sample}.sorted.bam" \
            > "$REGION_OUTPUT/bams/${sample}.flagstat.txt"

    done < "$REGION_OUTPUT/temp/diploid_samples.txt"
    
    cat > "$REGION_OUTPUT/region_summary.txt" <<EOF
Region: ${TARGET_CHR}/${REGION_NAME}
Date: $(date)
Alleles: $ALLELES_FASTA
Diploid samples: $NDIPLOID
BAM files: $REGION_OUTPUT/bams/
EOF
    
    echo ""
    
done < "$REGION_MAP"

sort -u "$ALL_SAMPLES_FILE" > "$UNIQUE_SAMPLES_FILE"
UNIQUE_SAMPLE_COUNT=$(wc -l < "$UNIQUE_SAMPLES_FILE" | awk '{print $1}')

if [[ -n "$CONTROL_REGION" ]]; then
    echo "Processing control region: ${CONTROL_REGION}"
    CONTROL_OUTPUT="$OUTPUT_DIR/${CONTROL_CHR}/${CONTROL_NAME}"
    mkdir -p "$CONTROL_OUTPUT"/{logs,temp,bams}

    CHROM_REFERENCE="$OUTPUT_DIR/${CONTROL_CHR}/reference_${CONTROL_CHR}.fa"
    if ! have_extracted_chr "$CONTROL_CHR"; then
        echo "Extracting chromosome ${CONTROL_CHR}..."
        mkdir -p "$OUTPUT_DIR/${CONTROL_CHR}"
        samtools faidx "$REFERENCE" "$CONTROL_CHR" > "$CHROM_REFERENCE"

        if [[ ! -s "$CHROM_REFERENCE" ]]; then
            echo "ERROR: Failed to extract $CONTROL_CHR, skipping control region"
            CONTROL_REGION=""
        else
            echo "Indexing ${CONTROL_CHR} with BWA..."
            bwa index "$CHROM_REFERENCE"
            samtools faidx "$CHROM_REFERENCE"
            mark_extracted_chr "$CONTROL_CHR"
        fi
    fi

    if [[ -n "$CONTROL_REGION" ]]; then
        CONTROL_FASTA="$CONTROL_OUTPUT/temp/${CONTROL_NAME}.reference.pansn.fa"
        CONTROL_GARGAMMEL_FASTA="$CONTROL_OUTPUT/temp/${CONTROL_NAME}.reference.fa"
        samtools faidx "$REFERENCE" "$CONTROL_REGION" > "$CONTROL_FASTA"
        sed "s/^>.*/>${CONTROL_NAME}#1#${CONTROL_CONTIG}/" "$CONTROL_FASTA" > "$CONTROL_OUTPUT/${CONTROL_NAME}.pansn.fa"
        sed "s/^>.*/>chr_endo/" "$CONTROL_FASTA" > "$CONTROL_GARGAMMEL_FASTA"

        cat >> "$GLOBAL_SUMMARY" <<EOF
${CONTROL_CHR}/${CONTROL_NAME}: haploid reference control
  Region: $CONTROL_REGION
  Output: $CONTROL_OUTPUT/bams/

EOF

        while IFS= read -r sample; do
            [[ -z "$sample" ]] && continue
            echo "  Control sample: ${sample}"
            SAMPLE_DIR="$CONTROL_OUTPUT/temp/${sample}"
            mkdir -p "$SAMPLE_DIR"/{endo,cont,bact}

            MAPPING_LOG="$CONTROL_OUTPUT/logs/${sample}_sequence_mapping.txt"
            echo "Sample: ${sample} | Control region: ${CONTROL_REGION}" > "$MAPPING_LOG"
            echo "Date: $(date)" >> "$MAPPING_LOG"
            echo "" >> "$MAPPING_LOG"

            cp "$CONTROL_GARGAMMEL_FASTA" "$SAMPLE_DIR/endo/control.fa"
            echo "[ENDOGENOUS CONTROL]" >> "$MAPPING_LOG"
            echo "control.fa <- $CONTROL_REGION from $REFERENCE" >> "$MAPPING_LOG"
            echo "haploid reference control; gargammel -c coverage is computed from this region" >> "$MAPPING_LOG"
            echo "" >> "$MAPPING_LOG"

            if (( $(echo "$CONT_RATIO > 0" | bc -l) )); then
                sed "s/^>.*/>chr_cont/" "$CONTROL_FASTA" > "$SAMPLE_DIR/cont/control_cont.fa"
                echo "[CONTAMINANT CONTROL] - reference control" >> "$MAPPING_LOG"
                echo "control_cont.fa <- $CONTROL_REGION from $REFERENCE" >> "$MAPPING_LOG"
                echo "" >> "$MAPPING_LOG"
            else
                sed "s/^>.*/>chr_cont/" "$CONTROL_FASTA" > "$SAMPLE_DIR/cont/dummy.fa"
                echo "[CONTAMINANT] - None" >> "$MAPPING_LOG"
                echo "" >> "$MAPPING_LOG"
            fi

            sed "s/^>.*/>chr_bact/" "$CONTROL_FASTA" > "$SAMPLE_DIR/bact/dummy.fa"
            echo "[BACTERIAL] - None" >> "$MAPPING_LOG"
            echo "" >> "$MAPPING_LOG"

            GARGAMMEL_RUN="$GARGAMMEL_CMD \
                -c $COVERAGE \
                --comp 0,$CONT_RATIO,$ENDO_RATIO \
                -l $FRAGMENT_LENGTH \
                -rl $READ_LENGTH \
                -o $SAMPLE_DIR/sim"

            if [[ -n "$DEAM_RATE" ]]; then
                GARGAMMEL_RUN="$GARGAMMEL_RUN -damage $DEAM_RATE"
            else
                if [[ "$DEAMINATION" == "single" ]]; then
                    MATRIX="$MATRIX_SINGLE"
                else
                    MATRIX="$MATRIX_DOUBLE"
                fi
                GARGAMMEL_RUN="$GARGAMMEL_RUN -matfile $MATRIX"
            fi

            if [[ "$LIBRARY_TYPE" == "se" ]]; then
                GARGAMMEL_RUN="$GARGAMMEL_RUN -se"
            fi

            GARGAMMEL_RUN="$GARGAMMEL_RUN $SAMPLE_DIR"
            eval $GARGAMMEL_RUN > "$CONTROL_OUTPUT/logs/${sample}_gargammel.log" 2>&1

            if [[ "$LIBRARY_TYPE" == "pe" ]]; then
                fastp \
                    --in1 "$SAMPLE_DIR/sim_s1.fq.gz" \
                    --in2 "$SAMPLE_DIR/sim_s2.fq.gz" \
                    --out1 "$SAMPLE_DIR/sim_s1_trimmed.fq.gz" \
                    --out2 "$SAMPLE_DIR/sim_s2_trimmed.fq.gz" \
                    --detect_adapter_for_pe \
                    --thread "$THREADS" \
                    --length_required 25 \
                    --json "$CONTROL_OUTPUT/logs/${sample}_fastp.json" \
                    --html "$CONTROL_OUTPUT/logs/${sample}_fastp.html" \
                    2> "$CONTROL_OUTPUT/logs/${sample}_fastp.log"

                READ1="$SAMPLE_DIR/sim_s1_trimmed.fq.gz"
                READ2="$SAMPLE_DIR/sim_s2_trimmed.fq.gz"
            else
                fastp \
                    --in1 "$SAMPLE_DIR/sim_s.fq.gz" \
                    --out1 "$SAMPLE_DIR/sim_s_trimmed.fq.gz" \
                    --thread "$THREADS" \
                    --length_required 25 \
                    --json "$CONTROL_OUTPUT/logs/${sample}_fastp.json" \
                    --html "$CONTROL_OUTPUT/logs/${sample}_fastp.html" \
                    --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAG \
                    2> "$CONTROL_OUTPUT/logs/${sample}_fastp.log"

                READ1="$SAMPLE_DIR/sim_s_trimmed.fq.gz"
            fi

            if [[ "$LIBRARY_TYPE" == "pe" ]]; then
                bwa aln -l 16500 -n 0.01 -o 2 -t "$THREADS" \
                    "$CHROM_REFERENCE" "$READ1" \
                    > "$SAMPLE_DIR/sim_1.sai" \
                    2> "$CONTROL_OUTPUT/logs/${sample}_bwa_aln_1.log"

                bwa aln -l 16500 -n 0.01 -o 2 -t "$THREADS" \
                    "$CHROM_REFERENCE" "$READ2" \
                    > "$SAMPLE_DIR/sim_2.sai" \
                    2> "$CONTROL_OUTPUT/logs/${sample}_bwa_aln_2.log"

                bwa sampe "$CHROM_REFERENCE" \
                    "$SAMPLE_DIR/sim_1.sai" "$SAMPLE_DIR/sim_2.sai" \
                    "$READ1" "$READ2" \
                    2> "$CONTROL_OUTPUT/logs/${sample}_bwa_sampe.log" | \
                samtools sort -@ "$THREADS" \
                    -o "$CONTROL_OUTPUT/bams/${sample}.sorted.bam" -
            else
                bwa aln -l 16500 -n 0.01 -o 2 -t "$THREADS" \
                    "$CHROM_REFERENCE" "$READ1" \
                    > "$SAMPLE_DIR/sim.sai" \
                    2> "$CONTROL_OUTPUT/logs/${sample}_bwa_aln.log"

                bwa samse "$CHROM_REFERENCE" \
                    "$SAMPLE_DIR/sim.sai" "$READ1" \
                    2> "$CONTROL_OUTPUT/logs/${sample}_bwa_samse.log" | \
                samtools sort -@ "$THREADS" \
                    -o "$CONTROL_OUTPUT/bams/${sample}.sorted.bam" -
            fi

            samtools index "$CONTROL_OUTPUT/bams/${sample}.sorted.bam"
            samtools flagstat "$CONTROL_OUTPUT/bams/${sample}.sorted.bam" \
                > "$CONTROL_OUTPUT/bams/${sample}.flagstat.txt"
        done < "$UNIQUE_SAMPLES_FILE"

        cat > "$CONTROL_OUTPUT/region_summary.txt" <<EOF
Region: ${CONTROL_CHR}/${CONTROL_NAME}
Date: $(date)
Reference control region: $CONTROL_REGION
Haploid control FASTA: $CONTROL_OUTPUT/${CONTROL_NAME}.pansn.fa
BAM files: $CONTROL_OUTPUT/bams/
EOF
    fi

    echo ""
fi

echo "" >> "$MERGE_SCRIPT"

while IFS= read -r sample; do
    [[ -z "$sample" ]] && continue
    cat >> "$MERGE_SCRIPT" <<EOF

echo "Merging: $sample"
SAMPLE_BAMS=()
EOF
    
    while IFS=$'\t' read -r chr region _; do
        [[ "$chr" =~ ^#.*$ ]] && continue
        [[ -z "$chr" || -z "$region" ]] && continue
        
        cat >> "$MERGE_SCRIPT" <<EOF
[[ -f "\$OUTPUT_DIR/${chr}/${region}/bams/${sample}.sorted.bam" ]] && SAMPLE_BAMS+=("\$OUTPUT_DIR/${chr}/${region}/bams/${sample}.sorted.bam")
EOF
    done < "$REGION_MAP"

    if [[ -n "$CONTROL_REGION" ]]; then
        cat >> "$MERGE_SCRIPT" <<EOF
[[ -f "\$OUTPUT_DIR/${CONTROL_CHR}/${CONTROL_NAME}/bams/${sample}.sorted.bam" ]] && SAMPLE_BAMS+=("\$OUTPUT_DIR/${CONTROL_CHR}/${CONTROL_NAME}/bams/${sample}.sorted.bam")
EOF
    fi
    
    cat >> "$MERGE_SCRIPT" <<EOF

if [[ \${#SAMPLE_BAMS[@]} -gt 1 ]]; then
    samtools merge -@ 4 "\$MERGED_DIR/${sample}.merged.bam" "\${SAMPLE_BAMS[@]}"
    samtools index "\$MERGED_DIR/${sample}.merged.bam"
    echo "  \${#SAMPLE_BAMS[@]} BAMs -> ${sample}.merged.bam"
elif [[ \${#SAMPLE_BAMS[@]} -eq 1 ]]; then
    cp "\${SAMPLE_BAMS[0]}" "\$MERGED_DIR/${sample}.merged.bam"
    samtools index "\$MERGED_DIR/${sample}.merged.bam"
    echo "  1 BAM -> ${sample}.merged.bam"
else
    echo "  WARNING: No BAMs found"
fi
EOF
done < "$UNIQUE_SAMPLES_FILE"

cat >> "$MERGE_SCRIPT" <<'EOF'

echo ""
echo "=== Complete ==="
echo "Merged BAMs: $MERGED_DIR"
EOF

chmod +x "$MERGE_SCRIPT"

cat >> "$GLOBAL_SUMMARY" <<EOF

Summary:
=========
Regions processed: $REGION_NUM
Control region: ${CONTROL_REGION:-none}
Unique samples: $UNIQUE_SAMPLE_COUNT

To merge BAMs per sample:
  bash $MERGE_SCRIPT

EOF

echo "=== AncestralSim Complete ==="
echo "Regions: $REGION_NUM | Samples: $UNIQUE_SAMPLE_COUNT"
echo "Output: $OUTPUT_DIR"
echo "Merge: bash $MERGE_SCRIPT"
echo ""
