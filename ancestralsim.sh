#!/bin/bash
# ancestralsim: Ancient DNA simulation from pangenome haplotypes

# default
COVERAGE=0.5
FRAGMENT_LENGTH=45
READ_LENGTH=75
DEAMINATION="single"
DEAM_RATE=""
LIBRARY_TYPE="pe"
REFERENCE=""
OUTPUT_DIR="output"
THREADS=4
GARGAMMEL_DIR=""
CONT_RATIO=0.1  # default 10% exogenous contamination

usage() {
cat <<EOF
Usage: $0 -f <alleles.fasta> -r <reference.fa> -g <gargammel_dir> [options]

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

EOF
exit 1
}

# params parser
while [[ $# -gt 0 ]]; do
    case $1 in
        -f) ALLELES_FASTA="$2"; shift 2 ;;
        -r) REFERENCE="$2"; shift 2 ;;
        -g) GARGAMMEL_DIR="$2"; shift 2 ;;
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

# input validation
if [[ -z "$ALLELES_FASTA" || -z "$REFERENCE" || -z "$GARGAMMEL_DIR" ]]; then
    echo "Error: -f, -r, and -g are required."; usage
fi
if [[ ! -d "$GARGAMMEL_DIR" ]]; then
    echo "Error: Gargammel directory not found: $GARGAMMEL_DIR"; exit 1
fi
if [[ ! -f "$GARGAMMEL_DIR/gargammel.pl" ]]; then
    echo "Error: gargammel.pl not found in $GARGAMMEL_DIR"; exit 1
fi

# matrix verification
MATRIX_DIR="$GARGAMMEL_DIR/src/matrices"
MATRIX_SINGLE="$MATRIX_DIR/single-"
MATRIX_DOUBLE="$MATRIX_DIR/double-"
if [[ ! -d "$MATRIX_DIR" ]]; then
    echo "Error: Matrices directory missing: $MATRIX_DIR"; exit 1
fi

# library type
if [[ "$LIBRARY_TYPE" != "se" ]] && [[ "$LIBRARY_TYPE" != "pe" ]]; then
    echo "Error: Library type must be 'se' or 'pe'"; exit 1
fi

# deamination type
if [[ "$DEAMINATION" != "single" ]] && [[ "$DEAMINATION" != "double" ]]; then
    echo "Error: Deamination type must be 'single' or 'double'"; exit 1
fi

# contamination-ratio validity
if (( $(echo "$CONT_RATIO < 0" | bc -l) )) || (( $(echo "$CONT_RATIO > 1" | bc -l) )); then
    echo "Error: --cont-ratio must be between 0 and 1"; exit 1
fi

# endogenous (comp order: bact, cont, endo)
ENDO_RATIO=$(python3 -c "print(round(1.0 - $CONT_RATIO, 6))")

# dependencies
for cmd in samtools bwa; do
    if ! command -v $cmd &> /dev/null; then
        echo "Error: $cmd not found. Please install it."; exit 1
    fi
done

# gargamel executable in conda environment
if command -v gargammel &>/dev/null; then
    GARGAMMEL_CMD="gargammel"
else
    GARGAMMEL_CMD="$GARGAMMEL_DIR/gargammel.pl"
fi

mkdir -p "$OUTPUT_DIR"/{logs,temp,bams}

echo "=== AncestralSim: Ancient DNA Simulation Pipeline ==="
echo "Input alleles: $ALLELES_FASTA"
echo "Reference: $REFERENCE"
echo "Gargammel: $GARGAMMEL_CMD"
echo "Coverage: ${COVERAGE}x"
echo "Fragment length: ${FRAGMENT_LENGTH}bp"
echo "Read length: ${READ_LENGTH}bp"
echo "Library type: $LIBRARY_TYPE"
echo "Deamination: $DEAMINATION-stranded"
echo "Contamination: ${CONT_RATIO} (endogenous: ${ENDO_RATIO})"
echo ""

# indexing if missing
if [[ ! -f "${ALLELES_FASTA}.fai" ]]; then
    echo "Indexing alleles FASTA..."
    samtools faidx "$ALLELES_FASTA"
fi
FAI_FILE="${ALLELES_FASTA}.fai"

# get diploid samples
echo "Identifying diploid samples..."
awk -F'[#\t]' '{print $1}' "$FAI_FILE" | sort | uniq -c | \
awk '$1==2 {print $2}' > "$OUTPUT_DIR/temp/diploid_samples.txt"
NDIPLOID=$(wc -l < "$OUTPUT_DIR/temp/diploid_samples.txt")
echo "Found ${NDIPLOID} diploid samples"
echo ""

if [[ $NDIPLOID -eq 0 ]]; then
    echo "Error: No diploid samples found in input FASTA"; exit 1
fi

# further reference indexing
if [[ ! -f "${REFERENCE}.bwt" ]]; then
    echo "Indexing reference with BWA..."
    bwa index "$REFERENCE"
fi
if [[ ! -f "${REFERENCE}.fai" ]]; then
    echo "Indexing reference with samtools..."
    samtools faidx "$REFERENCE"
fi

# rename fasta headers for garfamel
rename_fasta_header() {
    local input_fasta=$1
    local output_fasta=$2
    local new_name=$3
    local original_name=$4

    if [[ "$input_fasta" == *.gz ]]; then
        zcat "$input_fasta" | sed "s/^>${original_name}/>${new_name}/" > "$output_fasta"
    else
        sed "s/^>${original_name}/>${new_name}/" "$input_fasta" > "$output_fasta"
    fi
}

# process samples
echo "[$(date)] Starting simulation for ${NDIPLOID} samples..."
echo ""

while IFS= read -r sample; do
    echo ">> Processing sample: ${sample}"
    SAMPLE_DIR="$OUTPUT_DIR/temp/${sample}"
    mkdir -p "$SAMPLE_DIR"/{endo,cont,bact}
    
    # mapping log file
    MAPPING_LOG="$OUTPUT_DIR/logs/${sample}_sequence_mapping.txt"
    echo "Sequence Name Mapping for Sample: ${sample}" > "$MAPPING_LOG"
    echo "Generated: $(date)" >> "$MAPPING_LOG"
    echo "" >> "$MAPPING_LOG"

    # get endogenous haplotypes
    echo "   Extracting endogenous haplotypes..."
    
    # original sequences
    ENDO_HAP1_ORIG=$(grep "^${sample}#1#" "$FAI_FILE" | cut -f1)
    ENDO_HAP2_ORIG=$(grep "^${sample}#2#" "$FAI_FILE" | cut -f1)
    
    # extract and rename
    samtools faidx "$ALLELES_FASTA" "$ENDO_HAP1_ORIG" | \
        sed "s/^>.*/>chr_endo/" > "$SAMPLE_DIR/endo/hap1.fa"
    samtools faidx "$ALLELES_FASTA" "$ENDO_HAP2_ORIG" | \
        sed "s/^>.*/>chr_endo/" > "$SAMPLE_DIR/endo/hap2.fa"
    
    # to log
    echo "[ENDOGENOUS]" >> "$MAPPING_LOG"
    echo "hap1.fa (chr_endo) <- $ENDO_HAP1_ORIG" >> "$MAPPING_LOG"
    echo "hap2.fa (chr_endo) <- $ENDO_HAP2_ORIG" >> "$MAPPING_LOG"
    echo "" >> "$MAPPING_LOG"

    # select and extract contaminants
    if (( $(echo "$CONT_RATIO > 0" | bc -l) )); then
        CONT_SAMPLE=$(grep -v "^${sample}$" "$OUTPUT_DIR/temp/diploid_samples.txt" | shuf -n 1)
        echo "   Extracting contaminant: ${CONT_SAMPLE} (ratio: ${CONT_RATIO})"
        
        # original contaminant sequence names
        CONT_HAP1_ORIG=$(grep "^${CONT_SAMPLE}#1#" "$FAI_FILE" | cut -f1)
        CONT_HAP2_ORIG=$(grep "^${CONT_SAMPLE}#2#" "$FAI_FILE" | cut -f1)
        
        # extract t and rename to common chromosome name
        samtools faidx "$ALLELES_FASTA" "$CONT_HAP1_ORIG" | \
            sed "s/^>.*/>chr_cont/" > "$SAMPLE_DIR/cont/cont_hap1.fa"
        samtools faidx "$ALLELES_FASTA" "$CONT_HAP2_ORIG" | \
            sed "s/^>.*/>chr_cont/" > "$SAMPLE_DIR/cont/cont_hap2.fa"
        
        # to log
        echo "[CONTAMINANT] - Sample: ${CONT_SAMPLE}" >> "$MAPPING_LOG"
        echo "cont_hap1.fa (chr_cont) <- $CONT_HAP1_ORIG" >> "$MAPPING_LOG"
        echo "cont_hap2.fa (chr_cont) <- $CONT_HAP2_ORIG" >> "$MAPPING_LOG"
        echo "" >> "$MAPPING_LOG"
    else
        echo "   No contamination (creating dummy cont file)"
        # first endogenous haplotype as dummy - not used
        samtools faidx "$ALLELES_FASTA" "$ENDO_HAP1_ORIG" | \
            sed "s/^>.*/>chr_cont/" > "$SAMPLE_DIR/cont/dummy.fa"
        
        echo "[CONTAMINANT] - None (dummy file)" >> "$MAPPING_LOG"
        echo "dummy.fa (chr_cont) <- $ENDO_HAP1_ORIG (dummy)" >> "$MAPPING_LOG"
        echo "" >> "$MAPPING_LOG"
    fi
    
    # bacterial with dummy haplotype
    samtools faidx "$ALLELES_FASTA" "$ENDO_HAP1_ORIG" | \
        sed "s/^>.*/>chr_bact/" > "$SAMPLE_DIR/bact/dummy.fa"
    
    echo "[BACTERIAL] - None (dummy file)" >> "$MAPPING_LOG"
    echo "dummy.fa (chr_bact) <- $ENDO_HAP1_ORIG (dummy)" >> "$MAPPING_LOG"
    echo "" >> "$MAPPING_LOG"

    # build gargamel command
    echo "   Running gargammel simulation..."
    
    GARGAMMEL_RUN="$GARGAMMEL_CMD \
        -c $COVERAGE \
        --comp 0,$CONT_RATIO,$ENDO_RATIO \
        -l $FRAGMENT_LENGTH \
        -rl $READ_LENGTH \
        -o $SAMPLE_DIR/sim"
    
    # deamination
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
    
    # library flag
    if [[ "$LIBRARY_TYPE" == "se" ]]; then
        GARGAMMEL_RUN="$GARGAMMEL_RUN -se"
    fi
    
    # input directory - should be last
    GARGAMMEL_RUN="$GARGAMMEL_RUN $SAMPLE_DIR"
    
    # execute
    eval $GARGAMMEL_RUN > "$OUTPUT_DIR/logs/${sample}_gargammel.log" 2>&1

    # alignment
    echo "   Aligning reads with BWA aln..."
    
    if [[ "$LIBRARY_TYPE" == "pe" ]]; then
        # paired
        bwa aln \
            -l 16500 \
            -n 0.01 \
            -o 2 \
            -t "$THREADS" \
            "$REFERENCE" \
            "$SAMPLE_DIR/sim_s1.fq.gz" \
            > "$SAMPLE_DIR/sim_1.sai" \
            2> "$OUTPUT_DIR/logs/${sample}_bwa_aln_1.log"
        
        bwa aln \
            -l 16500 \
            -n 0.01 \
            -o 2 \
            -t "$THREADS" \
            "$REFERENCE" \
            "$SAMPLE_DIR/sim_s2.fq.gz" \
            > "$SAMPLE_DIR/sim_2.sai" \
            2> "$OUTPUT_DIR/logs/${sample}_bwa_aln_2.log"
        
        echo "   Converting to BAM..."
        bwa sampe \
            "$REFERENCE" \
            "$SAMPLE_DIR/sim_1.sai" \
            "$SAMPLE_DIR/sim_2.sai" \
            "$SAMPLE_DIR/sim_s1.fq.gz" \
            "$SAMPLE_DIR/sim_s2.fq.gz" \
            2> "$OUTPUT_DIR/logs/${sample}_bwa_sampe.log" | \
        samtools sort -@ "$THREADS" -o "$OUTPUT_DIR/bams/${sample}.sorted.bam" -
    else
        # single
        bwa aln \
            -l 16500 \
            -n 0.01 \
            -o 2 \
            -t "$THREADS" \
            "$REFERENCE" \
            "$SAMPLE_DIR/sim_s.fq.gz" \
            > "$SAMPLE_DIR/sim.sai" \
            2> "$OUTPUT_DIR/logs/${sample}_bwa_aln.log"
        
        echo "   Converting to BAM..."
        bwa samse \
            "$REFERENCE" \
            "$SAMPLE_DIR/sim.sai" \
            "$SAMPLE_DIR/sim_s.fq.gz" \
            2> "$OUTPUT_DIR/logs/${sample}_bwa_samse.log" | \
        samtools sort -@ "$THREADS" -o "$OUTPUT_DIR/bams/${sample}.sorted.bam" -
    fi

    # index and stats
    samtools index "$OUTPUT_DIR/bams/${sample}.sorted.bam"
    samtools flagstat "$OUTPUT_DIR/bams/${sample}.sorted.bam" \
        > "$OUTPUT_DIR/bams/${sample}.flagstat.txt"
    
    echo "   Complete: $OUTPUT_DIR/bams/${sample}.sorted.bam"
    echo "   Sequence mapping saved: $MAPPING_LOG"
    echo ""

done < "$OUTPUT_DIR/temp/diploid_samples.txt"

# small report
echo "[$(date)] Generating summary report..."

cat > "$OUTPUT_DIR/simulation_summary.txt" <<EOF
AncestralSim Simulation Report
==============================
Date: $(date)
Input alleles: $ALLELES_FASTA
Reference: $REFERENCE
Gargammel directory: $GARGAMMEL_DIR

Simulation Parameters:
- Coverage: ${COVERAGE}x
- Fragment length: ${FRAGMENT_LENGTH}bp
- Read length: ${READ_LENGTH}bp
- Library type: $LIBRARY_TYPE
- Deamination type: $DEAMINATION-stranded
- Contamination ratio: ${CONT_RATIO} (endogenous: ${ENDO_RATIO})
- Number of diploid samples: $NDIPLOID

Output Files:
- BAM files: $OUTPUT_DIR/bams/*.sorted.bam
- Sequence mappings: $OUTPUT_DIR/logs/*_sequence_mapping.txt
- Gargammel logs: $OUTPUT_DIR/logs/*_gargammel.log

Output BAM files:
EOF

ls -lh "$OUTPUT_DIR/bams/"*.bam >> "$OUTPUT_DIR/simulation_summary.txt"

echo ""
echo "=== AncestralSim Complete ==="
echo "Output directory: $OUTPUT_DIR"
echo "BAM files: $OUTPUT_DIR/bams/"
echo "Sequence mappings: $OUTPUT_DIR/logs/*_sequence_mapping.txt"
echo "Summary: $OUTPUT_DIR/simulation_summary.txt"