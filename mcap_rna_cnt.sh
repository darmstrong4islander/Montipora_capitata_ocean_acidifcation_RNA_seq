#!/bin/bash

# SLURM job script header
#SBATCH -J mcap_rna_seq                               # Job name
#SBATCH -e mcap_rna_cnt.err                           # Error file
#SBATCH -N 1                                          # Number of nodes
#SBATCH -n 64                                         # Number of cores
#SBATCH -p bigmem                                     # Partition for high memory
#SBATCH --qos=highmem                                 # QOS for memory-intensive jobs
#SBATCH --mail-user=darmstrong4@islander.tamucc.edu   # Email notifications
#SBATCH --mail-type=BEGIN,END,FAIL                    # Notify on job begin, end, and failure
#SBATCH --time=96:00:00                               # Maximum runtime

module load fastqc fastp hisat2 samtools stringtie gffcompare

RAW_READS="./00_fastq"
FASTQC_DIR="./fastqc_reports"
CLEANED_DIR="./cleaned_reads"
ALIGN_DIR="./hisat2"
STRINGTIE_DIR="./stringtie"
REF_DIR="./reference_genome"
OUTPUT_DIR="./output"
ANNOTATION_FILE="$REF_DIR/Mcap.GFFannotation.putnam.gff"

# Build HISAT2 index
if [ ! -f "$REF_DIR/ref_index.1.ht2" ]; then
  echo "Building HISAT2 index..."
  hisat2-build "$REF_DIR/Mcap.genome_assembly.fasta" "$REF_DIR/ref_index" || { echo "Index building failed"; exit 1; }
else
  echo "HISAT2 index already exists."
fi

# Process each sample
for file in "$RAW_READS"/*_R1_001.fastq.gz; do
  BASE=$(basename "$file" _R1_001.fastq.gz)
  echo "Processing sample: $BASE"

  # FastQC on raw reads
  echo "Running FastQC for $BASE..."
  fastqc -o "$FASTQC_DIR" "$RAW_READS/${BASE}_R1_001.fastq.gz" "$RAW_READS/${BASE}_R2_001.fastq.gz" || exit 1

  # Quality trimming with Fastp
  echo "Running Fastp for $BASE..."
  fastp \
    --in1 "$RAW_READS/${BASE}_R1_001.fastq.gz" \
    --in2 "$RAW_READS/${BASE}_R2_001.fastq.gz" \
    --out1 "$CLEANED_DIR/${BASE}_R1_001.clean.fastq.gz" \
    --out2 "$CLEANED_DIR/${BASE}_R2_001.clean.fastq.gz" \
    --html "$CLEANED_DIR/${BASE}_fastp_report.html" || exit 1

  # Align cleaned reads with HISAT2
  echo "Running HISAT2 alignment for $BASE..."
  hisat2 -p 32 --rf --dta -x "$REF_DIR/ref_index" \
    -1 "$CLEANED_DIR/${BASE}_R1_001.clean.fastq.gz" \
    -2 "$CLEANED_DIR/${BASE}_R2_001.clean.fastq.gz" \
    -S "$ALIGN_DIR/${BASE}.sam" || exit 1

  # Convert SAM to sorted BAM
  echo "Converting SAM to BAM for $BASE..."
  samtools sort -@ 32 -o "$ALIGN_DIR/${BASE}.bam" "$ALIGN_DIR/${BASE}.sam" || exit 1
  rm "$ALIGN_DIR/${BASE}.sam"

  # Run StringTie for transcript assembly
  echo "Running StringTie for $BASE..."
  stringtie "$ALIGN_DIR/${BASE}.bam" \
    -G "$ANNOTATION_FILE" \
    -o "$STRINGTIE_DIR/${BASE}.gtf" \
    -A "$STRINGTIE_DIR/${BASE}_gene_abundance.tab" \
    -p 32 --rf -e || exit 1

  # Remove BAM file to save space
  rm "$ALIGN_DIR/${BASE}.bam"

done
