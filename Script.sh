#!/usr/bin/env bash
###############################################################################
# DESeq2 Analysis Pipeline - Bash Preprocessing
#
# Consolidates the workflow described in:
#   1_Data_Collection.md   -> download raw fastq files
#   2_Quality_Control.md   -> FastQC + MultiQC
#   3_Trimming.md          -> Trim Galore
#   4_Alignment.md         -> HISAT2 alignment + samtools sort/index
#   5_Feature_Count.md     -> featureCounts count matrix
#
# Usage:
#   ./pipeline.sh all                # run every step in order
#   ./pipeline.sh download           # run only one step
#   ./pipeline.sh qc trim align      # run a specific subset, in order given
#
# Edit the CONFIG block below to match your project before running.
###############################################################################

set -euo pipefail

# ============================== CONFIG ======================================
THREADS=4

RAW_DIR="raw_fastq"                       # where downloaded fastq.gz files live
QC_DIR="fastqc_results"                   # FastQC / MultiQC output
TRIM_DIR="trimmed_fastq"                  # Trim Galore output
ALIGN_DIR="aligned"                       # BAM output
COUNTS_DIR="counts"                       # featureCounts output

DOWNLOAD_SCRIPT="download_links.sh"       # file of wget links, see 1_Data_Collection.md
GENOME_INDEX="genome_index/GRCh38_index"  # prefix passed to `hisat2 -x`
GTF_FILE="annotation.gtf"                 # reference annotation for featureCounts

TRIM_QUALITY=20
TRIM_LENGTH=30
# =============================================================================

log() { echo -e "\n>>> $*\n"; }

require() {
    command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found on PATH." >&2; exit 1; }
}

setup_dirs() {
    mkdir -p "$RAW_DIR" "$QC_DIR" "$TRIM_DIR" "$ALIGN_DIR" "$COUNTS_DIR"
}

# ---- Step 1: Data Collection (1_Data_Collection.md) ------------------------
step_download() {
    require wget
    log "Step 1: Downloading raw files via $DOWNLOAD_SCRIPT"
    [ -f "$DOWNLOAD_SCRIPT" ] || { echo "ERROR: $DOWNLOAD_SCRIPT not found." >&2; exit 1; }
    chmod +x "$DOWNLOAD_SCRIPT"
    ( cd "$RAW_DIR" && "$OLDPWD/$DOWNLOAD_SCRIPT" )
}

# ---- Step 2: Quality Control (2_Quality_Control.md) -------------------------
step_qc() {
    require fastqc
    require multiqc
    log "Step 2: Running FastQC on raw reads"
    fastqc "$RAW_DIR"/*.fastq.gz -o "$QC_DIR" -t "$THREADS"
    log "Step 2: Aggregating with MultiQC"
    multiqc "$QC_DIR" --interactive -o "$QC_DIR"
}

# ---- Step 3: Trimming (3_Trimming.md) ---------------------------------------
step_trim() {
    require trim_galore
    log "Step 3: Trimming paired-end reads with Trim Galore"
    shopt -s nullglob
    for file1 in "$RAW_DIR"/*_R1.fastq.gz; do
        base=$(basename "$file1" "_R1.fastq.gz")
        file2="$RAW_DIR/${base}_R2.fastq.gz"
        if [ ! -f "$file2" ]; then
            echo "WARNING: mate file $file2 not found, skipping $base" >&2
            continue
        fi
        echo "Trimming $base ..."
        trim_galore --paired --quality "$TRIM_QUALITY" --length "$TRIM_LENGTH" \
            --fastqc --output_dir "$TRIM_DIR" "$file1" "$file2"
    done
    shopt -u nullglob
}

# ---- Step 4: Alignment (4_Alignment.md) -------------------------------------
step_align() {
    require hisat2
    require samtools
    log "Step 4: Aligning trimmed reads with HISAT2"
    shopt -s nullglob
    for r1 in "$TRIM_DIR"/*_R1_val_1.fq.gz; do
        base=$(basename "$r1" "_R1_val_1.fq.gz")
        r2="$TRIM_DIR/${base}_R2_val_2.fq.gz"
        if [ ! -f "$r2" ]; then
            echo "WARNING: mate file $r2 not found, skipping $base" >&2
            continue
        fi

        echo "Aligning $base ..."
        hisat2 -p "$THREADS" -x "$GENOME_INDEX" -1 "$r1" -2 "$r2" \
            -S "$ALIGN_DIR/${base}.sam"

        echo "Converting/sorting/indexing $base ..."
        samtools view -bS "$ALIGN_DIR/${base}.sam" \
            | samtools sort -@ "$THREADS" -o "$ALIGN_DIR/${base}_sorted.bam"
        samtools index "$ALIGN_DIR/${base}_sorted.bam"
        rm -f "$ALIGN_DIR/${base}.sam"
    done
    shopt -u nullglob
}

# ---- Step 5: Feature Counts (5_Feature_Count.md) ----------------------------
step_featurecounts() {
    require featureCounts
    log "Step 5: Generating count matrix with featureCounts"
    featureCounts -T "$THREADS" -p -a "$GTF_FILE" \
        -o "$COUNTS_DIR/counts.txt" "$ALIGN_DIR"/*_sorted.bam
    log "Done. Count matrix written to $COUNTS_DIR/counts.txt"
    echo "Feed this file (and your sample metadata) into DGE_analysis.R"
}

# ---- Driver ------------------------------------------------------------------
run_step() {
    case "$1" in
        download)      step_download ;;
        qc)             step_qc ;;
        trim)           step_trim ;;
        align)          step_align ;;
        featurecounts)  step_featurecounts ;;
        all)
            step_download
            step_qc
            step_trim
            step_align
            step_featurecounts
            ;;
        *)
            echo "Unknown step: $1" >&2
            echo "Valid steps: download qc trim align featurecounts all" >&2
            exit 1
            ;;
    esac
}

main() {
    setup_dirs
    if [ "$#" -eq 0 ]; then
        echo "Usage: $0 [all|download|qc|trim|align|featurecounts] ..." >&2
        exit 1
    fi
    for step in "$@"; do
        run_step "$step"
    done
}

main "$@"
