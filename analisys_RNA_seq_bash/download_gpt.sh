#!/bin/bash


#acive with 
# chmod +x process_bam.sh
#./process_bam.sh
#
# Usage: ./rnaseq_pipeline.sh /path/to/fastq /path/to/genome.fa /path/to/annotation.gtf /output/dir

# Arguments
FASTQ_DIR=$1
GENOME_FA=$2
GTF_FILE=$3
OUT_DIR=$4
THREADS=8

# STAR genome index directory
GENOME_DIR="${OUT_DIR}/star_index"

mkdir -p $OUT_DIR $GENOME_DIR "${OUT_DIR}/qc" "${OUT_DIR}/alignments" "${OUT_DIR}/counts"

# Step 1: Index Genome with STAR (if not done)
if [ ! -f "${GENOME_DIR}/SA" ]; then
  echo "Indexing genome with STAR..."
  STAR --runThreadN $THREADS \
       --runMode genomeGenerate \
       --genomeDir $GENOME_DIR \
       --genomeFastaFiles $GENOME_FA \
       --sjdbGTFfile $GTF_FILE \
       --sjdbOverhang 100
fi

# Step 2: Loop through FASTQ files
for R1 in ${FASTQ_DIR}/*_R1.fastq.gz; do
  SAMPLE=$(basename $R1 _R1.fastq.gz)
  R2=${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz

  echo "Processing sample: $SAMPLE"

  # Step 3: Quality Check
  fastqc -t $THREADS -o "${OUT_DIR}/qc" $R1 $R2

  # Step 4: Trimming (optional)
  # fastp -i $R1 -I $R2 -o ${OUT_DIR}/${SAMPLE}_R1.trimmed.fq.gz -O ${OUT_DIR}/${SAMPLE}_R2.trimmed.fq.gz

  # Step 5: Alignment
  STAR --runThreadN $THREADS \
       --genomeDir $GENOME_DIR \
       --readFilesIn $R1 $R2 \
       --readFilesCommand zcat \
       --outFileNamePrefix ${OUT_DIR}/alignments/${SAMPLE}_ \
       --outSAMtype BAM SortedByCoordinate

  # Step 6: Index BAM
  samtools index ${OUT_DIR}/alignments/${SAMPLE}_Aligned.sortedByCoord.out.bam

done

# Step 7: Quantification with featureCounts
BAMS=$(ls ${OUT_DIR}/alignments/*_Aligned.sortedByCoord.out.bam | tr '\n' ' ')
featureCounts -T $THREADS -a $GTF_FILE -o ${OUT_DIR}/counts/gene_counts.txt $BAMS

echo "RNA-seq pipeline completed successfully!"
