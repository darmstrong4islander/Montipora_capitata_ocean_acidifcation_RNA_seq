#!/bin/bash

# SLURM job script header

#SBATCH -J mcap_data_cleaning                        # Job name
#SBATCH -e mcap_error_file.err                       # File for job errors
#SBATCH -N 1                                         # Number of nodes
#SBATCH -n 64                                        # Number of cores
#SBATCH -p bigmem                                    # Partition to use
#SBATCH --qos=highmem                                # QOS for memory-intensive jobs
#SBATCH --mail-user=darmstrong4@islander.tamucc.edu  # Email notifications
#SBATCH --mail-type=BEGIN,END,FAIL                   # Notify on job begin, end, and failure
#SBATCH --time=96:00:00                              # Maximum runtime

# Paths and environment variables
STRINGTIE_DIR="./stringtie"
OUTPUT_DIR="./output"
ANNOTATION_FILE="./reference_genome/Mcap.GFFannotation.putnam.gff"

# Step 1: Generate the list of GTF files for StringTie merge (file paths only)
MERGE_LIST="$STRINGTIE_DIR/merge_list.txt"
> "$MERGE_LIST"  # Clear the file before writing

for gtf_file in "$STRINGTIE_DIR"/*.gtf; do
  # Exclude merged GTF files
  if [[ $(basename "$gtf_file") != "merged.annotated.gtf" && $(basename "$gtf_file") != "stringtie_merged.gtf" ]]; then
    echo "$gtf_file" >> "$MERGE_LIST"
  fi
done

# Verify the generated merge_list.txt
echo "Generated merge_list.txt for StringTie merge:"
cat "$MERGE_LIST"

# Step 2: Run StringTie merge
echo "Merging GTF files with StringTie..."
stringtie --merge -p 8 -G "$ANNOTATION_FILE" -o "$STRINGTIE_DIR/stringtie_merged.gtf" "$MERGE_LIST" 2>> mcap_error_file.err
if [ $? -ne 0 ]; then
  echo "Error during StringTie merge. Check mcap_error_file.err for details."
  exit 1
fi

# Step 3: Assess assembly performance with gffcompare
echo "Assessing assembly performance with gffcompare..."
gffcompare -r "$ANNOTATION_FILE" -G -o "$STRINGTIE_DIR/merged" "$STRINGTIE_DIR/stringtie_merged.gtf" 2>> mcap_error_file.err
if [ $? -ne 0 ]; then
  echo "Error during GFF comparison. Check mcap_error_file.err for details."
  exit 1
fi

# Step 4: Generate the list of GTF files for prepDE.py3 (with sample IDs and paths)
PREPDE_LIST="$STRINGTIE_DIR/mergelist.txt"
> "$PREPDE_LIST"  # Clear the file before writing

for gtf_file in "$STRINGTIE_DIR"/*.gtf; do
  # Exclude merged GTF files
  if [[ $(basename "$gtf_file") != "merged.annotated.gtf" && $(basename "$gtf_file") != "stringtie_merged.gtf" ]]; then
    SAMPLE_NAME=$(basename "$gtf_file" .gtf)
    echo -e "${SAMPLE_NAME}\t${gtf_file}" >> "$PREPDE_LIST"
  fi
done

# Verify the generated mergelist.txt for prepDE.py3
echo "Generated mergelist.txt for prepDE.py3:"
cat "$PREPDE_LIST"

# Step 5: Generate Gene and Transcript Count Matrices
echo "Generating gene and transcript count matrices..."
./stringtie-2.2.1.Linux_x86_64/prepDE.py3 -i "$PREPDE_LIST" -g "$OUTPUT_DIR/gene_count_matrix.csv" -t "$OUTPUT_DIR/transcript_count_matrix.csv" 2>> mcap_error_file.err
if [ $? -ne 0 ]; then
  echo "Error during matrix generation. Check mcap_error_file.err for details."
  exit 1
fi

# End of script
echo "RNA-seq pipeline completed successfully."
