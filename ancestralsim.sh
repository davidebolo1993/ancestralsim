#!/bin/bash
# AncestralSim: Ancient DNA simulation from pangenome haplotypes
# This script simulates ancient DNA reads from pangenome graph alleles
# following best practices for aDNA simulation and alignment

set -e
set -o pipefail

# Default parameters following ancient DNA best practices
COVERAGE=0.5          # Low coverage typical of aDNA
FRAGMENT_LENGTH=45    # Mean fragment length (bp) - typical for aDNA [1]
FRAGMENT_SCALE=""     # Optional: scale parameter for log-normal distribution
FRAGMENT_LOC=""       # Optional: location parameter for log-normal distribution
READ_LENGTH=75        # Illumina read length
DEAMINATION="single"  # single or double strand deamination type
DEAM_RATE=""          # Optional: custom deamination rate (e.g., "0.03,0.4,0.01,0.3")
LIBRARY_TYPE="pe"     # pe (paired-end) or se (single-end)
REFERENCE=""          # Reference chromosome for alignment
OUTPUT_DIR="output"
THREADS=4
GARGAMMEL_DIR=""      # Path to gargammel installation directory

usage() {
    cat <<EOF
Usage: $0 -f <alleles.fasta> -r <reference.fa> -g <gargammel_dir> [options]

Required arguments:
  -f FILE    Input FASTA file with pangenome alleles (can be gzipped)
  -r FILE    Reference chromosome FASTA for alignment
  -g DIR     Path to gargammel installation directory

Optional arguments:
  -c FLOAT   Coverage (default: 0.5)
  -l INT     Mean fragment length in bp (default: 45)
             Typical ancient DNA: 30-100bp. Use with --loc/--scale for distribution
  --loc FLOAT   Location parameter for log-normal fragment distribution
  --scale FLOAT Scale parameter for log-normal fragment distribution
                Example: --loc 4.1 --scale 0.36 for realistic aDNA
  -d TYPE    Deamination type: single|double (default: single)
             single = single-stranded damage pattern (both ends damaged)
             double = double-stranded damage pattern (terminal damage only)
  --deam-rate "VALS"  Custom deamination rates (format: "5p_rate,5p_prob,3p_rate,3p_prob")
                      Example: "0.03,0.4,0.01,0.3" for 3% 5' and 1% 3' deamination
                      Higher values = older/more damaged samples
  -L TYPE    Library type: se|pe (default: pe)
  -R INT     Read length (default: 75)
  -o DIR     Output directory (default: output)
  -t INT     Number of threads (default: 4)
  -h         Show this help message

Fragment Size Notes:
  Ancient DNA is highly fragmented (typically 30-100bp). Fragmentation happens rapidly
  after death and stabilizes. Use -l for simple mean, or --loc/--scale for log-normal
  distribution matching real ancient samples [1,2].

Deamination Notes:
  C-to-T deamination at fragment ends is the hallmark of ancient DNA authenticity [3].
  Deamination increases with age and temperature. Use preset matrices (single/double)
  or specify custom rates with --deam-rate for specific scenarios [1,3].

References:
[1] Renaud et al. (2016) gargammel: a sequence simulator for ancient DNA. 
    Bioinformatics 33(4):577-579
[2] Sawyer et al. (2012) Temporal patterns of nucleotide misincorporations and DNA 
    fragmentation in ancient DNA. PLoS ONE 7(3):e34131
[3] Ginolhac et al. (2011) mapDamage: testing for damage patterns in ancient DNA 
    sequences. Bioinformatics 27(15):2153-2154
[4] Oliva et al. (2021) BWA-aln settings for ancient DNA alignment. 
    Ecol Evol 11(24):18743-18748
EOF
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -f) ALLELES_FASTA="$2"; shift 2 ;;
        -r) REFERENCE="$2"; shift 2 ;;
        -g) GARGAMMEL_DIR="$2"; shift 2 ;;
        -c) COVERAGE="$2"; shift 2 ;;
        -l) FRAGMENT_LENGTH="$2"; shift 2 ;;
        --loc) FRAGMENT_LOC="$2"; shift 2 ;;
        --scale) FRAGMENT_SCALE="$2"; shift 2 ;;
        -d) DEAMINATION="$2"; shift 2 ;;
        --deam-rate) DEAM_RATE="$2"; shift 2 ;;
        -L) LIBRARY_TYPE="$2"; shift 2 ;;
        -R) READ_LENGTH="$2"; shift 2 ;;
        -o) OUTPUT_DIR="$2"; shift 2 ;;
        -t) THREADS="$2"; shift 2 ;;
        -h) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

# Check required arguments
if [[ -z "$ALLELES_FASTA" ]] || [[ -z "$REFERENCE" ]] || [[ -z "$GARGAMMEL_DIR" ]]; then
    echo "Error: Missing required arguments"
    usage
fi

# Validate gargammel directory
if [[ ! -d "$GARGAMMEL_DIR" ]]; then
    echo "Error: Gargammel directory not found: $GARGAMMEL_DIR"
    exit 1
fi

# Check for gargammel.pl in the specified directory
if [[ ! -f "$GARGAMMEL_DIR/gargammel.pl" ]]; then
    echo "Error: gargammel.pl not found in $GARGAMMEL_DIR"
    echo "Please provide the correct gargammel installation directory"
    exit 1
fi

# Set matrix paths relative to gargammel installation
MATRIX_DIR="$GARGAMMEL_DIR/src/matrices"
if [[ ! -d "$MATRIX_DIR" ]]; then
    echo "Error: Matrices directory not found: $MATRIX_DIR"
    exit 1
fi

MATRIX_SINGLE="$MATRIX_DIR/single-"
MATRIX_DOUBLE="$MATRIX_DIR/double-"

# Verify that matrix files exist
if [[ -z "$DEAM_RATE" ]]; then
    if [[ "$DEAMINATION" == "single" ]]; then
        if ! ls ${MATRIX_SINGLE}* >/dev/null 2>&1; then
            echo "Error: Single-stranded deamination matrices not found at $MATRIX_SINGLE*"
            exit 1
        fi
    else
        if ! ls ${MATRIX_DOUBLE}* >/dev/null 2>&1; then
            echo "Error: Double-stranded deamination matrices not found at $MATRIX_DOUBLE*"
            exit 1
        fi
    fi
fi

# Validate library type
if [[ "$LIBRARY_TYPE" != "se" ]] && [[ "$LIBRARY_TYPE" != "pe" ]]; then
    echo "Error: Library type must be 'se' or 'pe'"
    exit 1
fi

# Validate deamination type
if [[ "$DEAMINATION" != "single" ]] && [[ "$DEAMINATION" != "double" ]]; then
    echo "Error: Deamination type must be 'single' or 'double'"
    exit 1
fi

# Check dependencies
for cmd in samtools bwa seqtk; do
    if ! command -v $cmd &> /dev/null; then
        echo "Error: $cmd not found. Please install it."
        exit 1
    fi
done

# Check for gargammel in PATH or use from specified directory
if command -v gargammel.pl &> /dev/null; then
    GARGAMMEL_CMD="gargammel.pl"
else
    GARGAMMEL_CMD="$GARGAMMEL_DIR/gargammel.pl"
fi

echo "=== AncestralSim: Ancient DNA Simulation Pipeline ==="
echo "Input alleles: $ALLELES_FASTA"
echo "Reference: $REFERENCE"
echo "Gargammel directory: $GARGAMMEL_DIR"
echo "Coverage: ${COVERAGE}x"
if [[ -n "$FRAGMENT_LOC" ]] && [[ -n "$FRAGMENT_SCALE" ]]; then
    echo "Fragment distribution: log-normal (loc=$FRAGMENT_LOC, scale=$FRAGMENT_SCALE)"
else
    echo "Fragment length: ${FRAGMENT_LENGTH}bp (mean)"
fi
echo "Read length: ${READ_LENGTH}bp"
echo "Library type: $LIBRARY_TYPE"
echo "Deamination: $DEAMINATION-stranded"
if [[ -n "$DEAM_RATE" ]]; then
    echo "Custom deamination rate: $DEAM_RATE"
else
    if [[ "$DEAMINATION" == "single" ]]; then
        echo "Deamination matrices: $MATRIX_SINGLE*"
    else
        echo "Deamination matrices: $MATRIX_DOUBLE*"
    fi
fi
echo ""

# Create output directory structure
mkdir -p "$OUTPUT_DIR"/{logs,temp,bams}

# Step 1: Extract diploid samples from pansn-formatted FASTA
echo "[$(date)] Step 1: Identifying diploid samples from pangenome alleles..."

# Parse FAI or create it if needed
FASTA_INPUT="$ALLELES_FASTA"
if [[ "$ALLELES_FASTA" == *.gz ]]; then
    if [[ ! -f "${ALLELES_FASTA}.fai" ]]; then
        samtools faidx "$ALLELES_FASTA"
    fi
    FAI_FILE="${ALLELES_FASTA}.fai"
else
    if [[ ! -f "${ALLELES_FASTA}.fai" ]]; then
        samtools faidx "$ALLELES_FASTA"
    fi
    FAI_FILE="${ALLELES_FASTA}.fai"
fi

# Extract sample names and identify diploid samples
# pansn format: sample#haplotype#contig_id
awk -F'[#\t]' '{print $1"#"$2}' "$FAI_FILE" | sort | uniq -c | \
    awk '$1==2 {print $2}' | cut -d'#' -f1 > "$OUTPUT_DIR/temp/diploid_samples.txt"

NDIPLOID=$(wc -l < "$OUTPUT_DIR/temp/diploid_samples.txt")
echo "Found $NDIPLOID diploid samples"

if [[ $NDIPLOID -eq 0 ]]; then
    echo "Error: No diploid samples found in input FASTA"
    exit 1
fi

# Index reference if needed
if [[ ! -f "${REFERENCE}.bwt" ]]; then
    echo "[$(date)] Indexing reference with BWA..."
    bwa index "$REFERENCE"
fi

if [[ ! -f "${REFERENCE}.fai" ]]; then
    echo "[$(date)] Indexing reference with samtools..."
    samtools faidx "$REFERENCE"
fi

# Step 2: Process each diploid sample
echo "[$(date)] Step 2: Simulating ancient DNA for each diploid sample..."

while IFS= read -r sample; do
    echo "  Processing sample: $sample"
    
    SAMPLE_DIR="$OUTPUT_DIR/temp/${sample}"
    mkdir -p "$SAMPLE_DIR"/{endo,cont,bact}
    
    # Extract both haplotypes for this sample
    # Haplotype 1
    if [[ "$FASTA_INPUT" == *.gz ]]; then
        samtools faidx "$FASTA_INPUT" $(grep "^${sample}#1#" "$FAI_FILE" | cut -f1) | \
            gzip -c > "$SAMPLE_DIR/endo/hap1.fa.gz"
    else
        samtools faidx "$FASTA_INPUT" $(grep "^${sample}#1#" "$FAI_FILE" | cut -f1) > \
            "$SAMPLE_DIR/endo/hap1.fa"
    fi
    
    # Haplotype 2
    if [[ "$FASTA_INPUT" == *.gz ]]; then
        samtools faidx "$FASTA_INPUT" $(grep "^${sample}#2#" "$FAI_FILE" | cut -f1) | \
            gzip -c > "$SAMPLE_DIR/endo/hap2.fa.gz"
    else
        samtools faidx "$FASTA_INPUT" $(grep "^${sample}#2#" "$FAI_FILE" | cut -f1) > \
            "$SAMPLE_DIR/endo/hap2.fa"
    fi
    
    # Create empty cont and bact directories (no contamination)
    touch "$SAMPLE_DIR/cont/.empty"
    touch "$SAMPLE_DIR/bact/.empty"
    
    # Step 3: Run gargammel simulation
    # Parameters follow best practices for aDNA simulation [1,2]:
    # - Fragment length: typically 30-100bp, often ~45bp mean [1,2]
    # - Deamination: C-to-T at 5' ends, G-to-A at 3' ends [1,3]
    # - Single-stranded vs double-stranded damage patterns [1]
    # - Coverage: typically <1x for ancient samples [1]
    echo "  Running gargammel simulation ($LIBRARY_TYPE, $DEAMINATION-stranded)..."
    
    # Build gargammel command
    GARGAMMEL_RUN="$GARGAMMEL_CMD \
        -c $COVERAGE \
        --comp 0,0,1 \
        -rl $READ_LENGTH \
        -o $SAMPLE_DIR/sim \
        $SAMPLE_DIR"
    
    # Add fragment length parameters
    if [[ -n "$FRAGMENT_LOC" ]] && [[ -n "$FRAGMENT_SCALE" ]]; then
        # Use log-normal distribution for realistic fragment size distribution
        GARGAMMEL_RUN="$GARGAMMEL_RUN --loc $FRAGMENT_LOC --scale $FRAGMENT_SCALE"
    else
        # Use simple mean fragment length
        GARGAMMEL_RUN="$GARGAMMEL_RUN -l $FRAGMENT_LENGTH"
    fi
    
    # Add deamination parameters
    if [[ -n "$DEAM_RATE" ]]; then
        # Use custom deamination rate
        GARGAMMEL_RUN="$GARGAMMEL_RUN -damage $DEAM_RATE"
    else
        # Use pre-calculated matrix
        if [[ "$DEAMINATION" == "single" ]]; then
            MATRIX="$MATRIX_SINGLE"
        else
            MATRIX="$MATRIX_DOUBLE"
        fi
        GARGAMMEL_RUN="$GARGAMMEL_RUN -matfile $MATRIX"
    fi
    
    # Add -se flag for single-end libraries
    if [[ "$LIBRARY_TYPE" == "se" ]]; then
        GARGAMMEL_RUN="$GARGAMMEL_RUN -se"
    fi
    
    eval $GARGAMMEL_RUN > "$OUTPUT_DIR/logs/${sample}_gargammel.log" 2>&1
    
    # Step 4: Align reads using BWA aln (best practice for aDNA) [4]
    # BWA aln is preferred over BWA-MEM for short, damaged aDNA reads
    # as it shows higher sensitivity for fragments <60bp with deamination [4]
    echo "  Aligning reads with BWA aln..."
    
    # BWA aln parameters optimized for ancient DNA [4]:
    # -l 16500: disable seed length threshold (sensitive to short reads)
    # -n 0.01: maximum edit distance
    # -o 2: gap open penalty
    
    if [[ "$LIBRARY_TYPE" == "pe" ]]; then
        # Paired-end alignment
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
        
        # Convert to BAM using bwa sampe
        echo "  Converting alignments to BAM..."
        bwa sampe \
            "$REFERENCE" \
            "$SAMPLE_DIR/sim_1.sai" \
            "$SAMPLE_DIR/sim_2.sai" \
            "$SAMPLE_DIR/sim_s1.fq.gz" \
            "$SAMPLE_DIR/sim_s2.fq.gz" \
            2> "$OUTPUT_DIR/logs/${sample}_bwa_sampe.log" | \
        samtools view -bS - | \
        samtools sort -@ "$THREADS" -o "$OUTPUT_DIR/bams/${sample}.sorted.bam" -
    else
        # Single-end alignment
        bwa aln \
            -l 16500 \
            -n 0.01 \
            -o 2 \
            -t "$THREADS" \
            "$REFERENCE" \
            "$SAMPLE_DIR/sim.fq.gz" \
            > "$SAMPLE_DIR/sim.sai" \
            2> "$OUTPUT_DIR/logs/${sample}_bwa_aln.log"
        
        # Convert to BAM using bwa samse
        echo "  Converting alignments to BAM..."
        bwa samse \
            "$REFERENCE" \
            "$SAMPLE_DIR/sim.sai" \
            "$SAMPLE_DIR/sim.fq.gz" \
            2> "$OUTPUT_DIR/logs/${sample}_bwa_samse.log" | \
        samtools view -bS - | \
        samtools sort -@ "$THREADS" -o "$OUTPUT_DIR/bams/${sample}.sorted.bam" -
    fi
    
    # Index the BAM file
    samtools index "$OUTPUT_DIR/bams/${sample}.sorted.bam"
    
    # Generate alignment statistics
    samtools flagstat "$OUTPUT_DIR/bams/${sample}.sorted.bam" \
        > "$OUTPUT_DIR/bams/${sample}.flagstat.txt"
    
    echo "  Sample $sample complete: $OUTPUT_DIR/bams/${sample}.sorted.bam"
    
done < "$OUTPUT_DIR/temp/diploid_samples.txt"

# Step 5: Generate summary report
echo "[$(date)] Step 5: Generating summary report..."

cat > "$OUTPUT_DIR/simulation_summary.txt" <<EOF
AncestralSim Simulation Report
==============================
Date: $(date)
Input alleles: $ALLELES_FASTA
Reference: $REFERENCE
Gargammel directory: $GARGAMMEL_DIR
Coverage: ${COVERAGE}x
EOF

if [[ -n "$FRAGMENT_LOC" ]] && [[ -n "$FRAGMENT_SCALE" ]]; then
    echo "Fragment distribution: log-normal (loc=$FRAGMENT_LOC, scale=$FRAGMENT_SCALE)" >> "$OUTPUT_DIR/simulation_summary.txt"
else
    echo "Fragment length: ${FRAGMENT_LENGTH}bp (mean)" >> "$OUTPUT_DIR/simulation_summary.txt"
fi

cat >> "$OUTPUT_DIR/simulation_summary.txt" <<EOF
Read length: ${READ_LENGTH}bp
Library type: $LIBRARY_TYPE
Deamination type: $DEAMINATION-stranded
EOF

if [[ -n "$DEAM_RATE" ]]; then
    echo "Custom deamination rate: $DEAM_RATE" >> "$OUTPUT_DIR/simulation_summary.txt"
else
    if [[ "$DEAMINATION" == "single" ]]; then
        echo "Deamination matrices: $MATRIX_SINGLE*" >> "$OUTPUT_DIR/simulation_summary.txt"
    else
        echo "Deamination matrices: $MATRIX_DOUBLE*" >> "$OUTPUT_DIR/simulation_summary.txt"
    fi
fi

cat >> "$OUTPUT_DIR/simulation_summary.txt" <<EOF
Number of diploid samples: $NDIPLOID

Output BAM files:
EOF

ls -lh "$OUTPUT_DIR/bams/"*.bam >> "$OUTPUT_DIR/simulation_summary.txt"

echo ""
echo "=== Simulation Complete ==="
echo "Output directory: $OUTPUT_DIR"
echo "BAM files: $OUTPUT_DIR/bams/"
echo "Summary: $OUTPUT_DIR/simulation_summary.txt"
echo ""
echo "References:"
echo "[1] Renaud et al. (2016) gargammel: a sequence simulator for ancient DNA"
echo "    Bioinformatics 33(4):577-579. doi:10.1093/bioinformatics/btw670"
echo "[2] Sawyer et al. (2012) Temporal patterns of nucleotide misincorporations"
echo "    and DNA fragmentation in ancient DNA. PLoS ONE 7(3):e34131"
echo "[3] Ginolhac et al. (2011) mapDamage: testing for damage patterns in aDNA"
echo "    Bioinformatics 27(15):2153-2154. doi:10.1093/bioinformatics/btr347"
echo "[4] Oliva et al. (2021) BWA-aln settings for ancient DNA alignment"
echo "    Ecol Evol 11(24):18743-18748. doi:10.1002/ece3.8297"

