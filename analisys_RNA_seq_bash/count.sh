# Download a test BAM file (human chr22 subset)
wget https://github.com/example-rnaseq/test-data/raw/main/sample.bam

# Download GTF (human chr22 from Ensembl)
wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.chr22.gtf.gz
gunzip Homo_sapiens.GRCh38.99.chr22.gtf.gz

htseq-count \
  --format bam \          # Input format (BAM/SAM)
  --order pos \          # BAM is position-sorted
  --stranded no \        # Library is unstranded (use 'yes' or 'reverse' if stranded)
  --type exon \          # Count reads overlapping exons
  --idattr gene_id \     # Use 'gene_id' in GTF as identifier
  sample.bam \           # Input BAM file
  Homo_sapiens.GRCh38.99.chr22.gtf > counts.txt