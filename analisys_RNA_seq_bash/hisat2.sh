#!/bin/bash

# HISAT2 Alignment Script
# Input: Paired-end FASTQ files
# Output: Aligned SAM/BAM files

# === Parameters ===
THREADS=4                         # Number of CPU threads
GENOME_INDEX="grch38/genome_snp_tran"  # HISAT2 index (built from reference)
FASTQ1="example_reads_1.fq"       # Forward reads
FASTQ2="example_reads_2.fq"       # Reverse reads
OUTPUT_SAM="output.sam"           # Output SAM file
OUTPUT_BAM="output.bam"           # Output BAM file (sorted)

# === Step 1: Run HISAT2 ===
hisat2 \
  -p $THREADS \                   # Threads for faster alignment
  -x $GENOME_INDEX \              # Reference genome index
  -1 $FASTQ1 \                    # Forward reads
  -2 $FASTQ2 \                    # Reverse reads
  -S $OUTPUT_SAM                  # Output SAM file

# === Step 2: Convert SAM to BAM & Sort ===
samtools view -@ $THREADS -b $OUTPUT_SAM | samtools sort -@ $THREADS -o $OUTPUT_BAM

# === Step 3: Index BAM (optional) ===
samtools index $OUTPUT_BAM

# === Cleanup ===
rm $OUTPUT_SAM  # Remove intermediate SAM file

echo "HISAT2 alignment complete. Output: $OUTPUT_BAM"