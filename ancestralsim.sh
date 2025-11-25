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

usage() {
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
  -o DIR      Output directory (default output)
  -t INT      Threads (default 4)
  -h          Show this help

EOF
exit 1
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
        -o) OUTPUT_DIR="$2"; shift 2 ;;
        -t) THREADS="$2"; shift 2 ;;
        -h) usage ;;
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
OUTPUT_DIR=$(realpath "$OUTPUT_DIR")
GARGAMMEL_DIR=$(realpath "$GARGAMMEL_DIR")
REGION_MAP=$(realpath "$REGION_MAP")

USE_SINGULARITY=false
SINGULARITY_CONTAINER=""

if command -v gargammel &>/dev/null; then
    GARGAMMEL_CMD="gargammel"
    echo "Using gargammel from PATH"
elif command -v singularity &>/dev/null; then
    SINGULARITY_CONTAINER="$GARGAMMEL_DIR/gargammel_1.1.4--hb66fcc3_0.sif"
    if [[ -f "$SINGULARITY_CONTAINER" ]]; then
        USE_SINGULARITY=true
        echo "Using gargammel via Singularity container"
        
        BIND_PATHS=""
        declare -A BIND_DIRS
        
        while IFS=$'\t' read -r _ _ alleles_file; do
            [[ "$alleles_file" =~ ^#.*$ ]] && continue
            [[ -z "$alleles_file" ]] && continue
            alleles_file=$(realpath "$alleles_file" 2>/dev/null) || continue
            ALLELES_DIR=$(dirname "$alleles_file")
            BIND_DIRS["$ALLELES_DIR"]=1
        done < "$REGION_MAP"
        
        REF_DIR=$(dirname "$REFERENCE")
        BIND_DIRS["$REF_DIR"]=1
        BIND_DIRS["$OUTPUT_DIR"]=1
        BIND_DIRS["$GARGAMMEL_DIR"]=1
        
        for dir in "${!BIND_DIRS[@]}"; do
            if [[ -z "$BIND_PATHS" ]]; then
                BIND_PATHS="$dir"
            else
                BIND_PATHS="$BIND_PATHS,$dir"
            fi
        done
        
        GARGAMMEL_CMD="singularity exec -B \"$BIND_PATHS\" $SINGULARITY_CONTAINER gargammel"
    else
        echo "Error: Singularity container not found: $SINGULARITY_CONTAINER"
        exit 1
    fi
else
    echo "Error: Neither gargammel nor singularity found."; exit 1
fi

if [[ "$USE_SINGULARITY" == false ]]; then
    if [[ ! -f "$GARGAMMEL_DIR/gargammel.pl" ]]; then
        echo "Error: gargammel.pl not found in $GARGAMMEL_DIR"; exit 1
    fi
fi

MATRIX_DIR="$GARGAMMEL_DIR/src/matrices"
MATRIX_SINGLE="$MATRIX_DIR/single-"
MATRIX_DOUBLE="$MATRIX_DIR/double-"

if [[ ! -d "$MATRIX_DIR" ]]; then
    echo "Error: Matrices directory missing: $MATRIX_DIR"; exit 1
fi
if [[ "$USE_SINGULARITY" == true ]]; then
    BIND_PATHS="$BIND_PATHS,$MATRIX_DIR"
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

for cmd in samtools bwa fastp; do
    if ! command -v $cmd &> /dev/null; then
        echo "Error: $cmd not found. Please install it."; exit 1
    fi
done

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
echo ""

declare -A CHR_EXTRACTED
N_REGIONS=$(grep -v "^#" "$REGION_MAP" | grep -v "^$" | wc -l)
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

declare -A ALL_SAMPLES

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
    
    if [[ -z "${CHR_EXTRACTED[$TARGET_CHR]}" ]]; then
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
        
        CHR_EXTRACTED[$TARGET_CHR]=1
    fi
    
    ALIGN_REFERENCE="$CHROM_REFERENCE"
    
    if [[ ! -f "${ALLELES_FASTA}.fai" ]]; then
        samtools faidx "$ALLELES_FASTA"
    fi
    FAI_FILE="${ALLELES_FASTA}.fai"
    
    # Identify diploid samples (2 alleles per sample)
    awk -F'[#\t]' '{print $1}' "$FAI_FILE" | sort | uniq -c | \
    awk '$1==2 {print $2}' > "$REGION_OUTPUT/temp/diploid_samples.txt"
    NDIPLOID=$(wc -l < "$REGION_OUTPUT/temp/diploid_samples.txt")
    
    if [[ $NDIPLOID -eq 0 ]]; then
        echo "WARNING: No diploid samples found, skipping"
        echo ""
        continue
    fi
    
    echo "Found ${NDIPLOID} diploid samples"
    echo ""
    
    while IFS= read -r sample; do
        ALL_SAMPLES[$sample]=1
    done < "$REGION_OUTPUT/temp/diploid_samples.txt"
    
    cat >> "$GLOBAL_SUMMARY" <<EOF
${TARGET_CHR}/${REGION_NAME}: $NDIPLOID samples
  Alleles: $ALLELES_FASTA
  Output: $REGION_OUTPUT/bams/

EOF
    
    while IFS= read -r sample; do
        echo "  Sample: ${sample}"
        SAMPLE_DIR="$REGION_OUTPUT/temp/${sample}"
        mkdir -p "$SAMPLE_DIR"/{endo,cont,bact}
        
        MAPPING_LOG="$REGION_OUTPUT/logs/${sample}_sequence_mapping.txt"
        echo "Sample: ${sample} | Chromosome: ${TARGET_CHR} | Region: ${REGION_NAME}" > "$MAPPING_LOG"
        echo "Date: $(date)" >> "$MAPPING_LOG"
        echo "" >> "$MAPPING_LOG"

        ENDO_HAP1_ORIG=$(grep "^${sample}#1#" "$FAI_FILE" | cut -f1)
        ENDO_HAP2_ORIG=$(grep "^${sample}#2#" "$FAI_FILE" | cut -f1)
        
        samtools faidx "$ALLELES_FASTA" "$ENDO_HAP1_ORIG" | \
            sed "s/^>.*/>chr_endo/" > "$SAMPLE_DIR/endo/hap1.fa"
        samtools faidx "$ALLELES_FASTA" "$ENDO_HAP2_ORIG" | \
            sed "s/^>.*/>chr_endo/" > "$SAMPLE_DIR/endo/hap2.fa"
        
        echo "[ENDOGENOUS]" >> "$MAPPING_LOG"
        echo "hap1.fa <- $ENDO_HAP1_ORIG" >> "$MAPPING_LOG"
        echo "hap2.fa <- $ENDO_HAP2_ORIG" >> "$MAPPING_LOG"
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

echo "" >> "$MERGE_SCRIPT"

for sample in "${!ALL_SAMPLES[@]}"; do
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
done

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
Unique samples: ${#ALL_SAMPLES[@]}

To merge BAMs per sample:
  bash $MERGE_SCRIPT

EOF

echo "=== AncestralSim Complete ==="
echo "Regions: $REGION_NUM | Samples: ${#ALL_SAMPLES[@]}"
echo "Output: $OUTPUT_DIR"
echo "Merge: bash $MERGE_SCRIPT"
echo ""

