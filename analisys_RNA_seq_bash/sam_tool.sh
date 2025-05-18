#!/bin/bash

# === Input BAM (replace with your file) ===
INPUT_BAM="example.bam"
SORTED_BAM="sorted.bam"
FILTERED_BAM="filtered.bam"
STATS_FILE="alignment_stats.txt"

# === Step 1: Sort BAM (required for indexing) ===
samtools sort -@ 4 -o $SORTED_BAM $INPUT_BAM

# === Step 2: Index BAM (for visualization/analysis) ===
samtools index $SORTED_BAM

# === Step 3: Filter BAM (optional) ===
# Remove unmapped reads (-F 4) and low-quality mappings (-q 20)
samtools view -@ 4 -b -F 4 -q 20 $SORTED_BAM > $FILTERED_BAM
samtools index $FILTERED_BAM

# === Step 4: Generate Alignment Stats ===
samtools flagstat $SORTED_BAM > $STATS_FILE
samtools idxstats $SORTED_BAM >> $STATS_FILE

echo "BAM processing complete. Outputs:"
echo "- Sorted BAM: $SORTED_BAM"
echo "- Filtered BAM: $FILTERED_BAM"
echo "- Stats: $STATS_FILE"