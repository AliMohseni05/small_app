┌───────────────────────────────────────────────────────┐
│                 RNA-seq Analysis Pipeline             │
└───────────────────────────────────────────────────────┘
                              │
                              ▼
┌───────────────────────────────────────────────────────┐
│ 1. Data Acquisition & Quality Control                 │
│    - Linux: SRA Toolkit (prefetch/fasterq-dump)       │
│    - FastQC (Quality Check)                           │
│    - Trimmomatic/Cutadapt (Adapter Trimming)          │
└───────────────────────────────────────────────────────┘
                              │
                              ▼
┌───────────────────────────────────────────────────────┐
│ 2. Alignment                                          │
│    - Linux: STAR/HISAT2 (Alignment to Reference)      │
│    - Samtools (BAM Processing)                        │
└───────────────────────────────────────────────────────┘
                              │
                              ▼
┌───────────────────────────────────────────────────────┐
│ 3. Transcript Assembly                                │
│    - Linux: StringTie/Cufflinks                       │
│    - GFFCompare (Transcript Comparison)               │
└───────────────────────────────────────────────────────┘
                              │
                              ▼
┌───────────────────────────────────────────────────────┐
│ 4. lncRNA Identification                              │
│    - CPC2/CNCI (Coding Potential)                     │
│    - FEELnc (lncRNA Classification)                   │
│    - Pfam Scan (Protein Domain Check)                 │
└───────────────────────────────────────────────────────┘
                              │
                              ▼
┌───────────────────────────────────────────────────────┐
│ 5. Quantification & DE Analysis                       │
│    - featureCounts (Read Counting)                    │
│    - R: DESeq2/edgeR (Differential Expression)        │
└───────────────────────────────────────────────────────┘
                              │
                              ▼
┌───────────────────────────────────────────────────────┐
│ 6. Functional Analysis                                 │
│    - R: clusterProfiler (GO/KEGG Enrichment)           │
│    - Python: GSEApy (Gene Set Enrichment)              │
│    - LncRNA2Target (Target Prediction)                │
└───────────────────────────────────────────────────────┘
                              │
                              ▼
┌───────────────────────────────────────────────────────┐
│ 7. Visualization & Reporting                          │
│    - R: ggplot2/ComplexHeatmap                        │
│    - Python: Matplotlib/Seaborn                       │
│    - IGV (Genome Browser)                             │
└───────────────────────────────────────────────────────┘
