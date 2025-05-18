#!/usr/bin/bash

# ----------
# Improved NGS Data Processing Pipeline
# ----------

# Exit on error and undefined variables for safety
set -euo pipefail

# ----- Configuration -----
# Use variables for flexibility
ACCESSION="SRR000001"
BASE_DIR="/home/skyland/my_atu"
ADAPTER_FILE="$BASE_DIR/adapters/TruSeq3-PE.fa"  # Ensure this file exists!

# Create directory structure
RAW_SRA="$BASE_DIR/sra_data"
RAW_FASTQ="$BASE_DIR/raw_fastq"
FASTQC_DIR="$BASE_DIR/fastqc_results"
TRIMMED_DIR="$BASE_DIR/trimmed_fastq"

mkdir -p "$RAW_SRA" "$RAW_FASTQ" "$FASTQC_DIR" "$TRIMMED_DIR"

# ----- Data Download -----
echo "Downloading $ACCESSION..."
prefetch "$ACCESSION" -O "$RAW_SRA"

# ----- SRA to FASTQ Conversion -----
echo "Converting to FASTQ..."
fasterq-dump "$ACCESSION" -O "$RAW_FASTQ" --split-files  # For paired-end data

# ----- Quality Control -----
echo "Running FastQC..."
fastqc "$RAW_FASTQ"/"$ACCESSION"* -o "$FASTQC_DIR"

# ----- Read Trimming -----
echo "Trimming reads..."

# Process paired-end files
R1="$RAW_FASTQ/${ACCESSION}_1.fastq"
R2="$RAW_FASTQ/${ACCESSION}_2.fastq"

trimmomatic PE -phred33 \
  "$R1" "$R2" \
  "$TRIMMED_DIR/${ACCESSION}_1_trimmed_paired.fastq" \
  "$TRIMMED_DIR/${ACCESSION}_1_trimmed_unpaired.fastq" \
  "$TRIMMED_DIR/${ACCESSION}_2_trimmed_paired.fastq" \
  "$TRIMMED_DIR/${ACCESSION}_2_trimmed_unpaired.fastq" \
  ILLUMINACLIP:"${ADAPTER_FILE}":2:30:10 \
  SLIDINGWINDOW:4:20 \
  MINLEN:25

echo "Processing complete!"