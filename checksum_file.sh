#!/bin/bash

# Script to verify that the files transferred correctly and generate a checksum report
# SLURM job script header (for a job submission)

#SBATCH -J mcap_checksum                             # Name of the job
#SBATCH -e checksum_file.err                         # File for job errors
#SBATCH -N 1                                         # Number of nodes
#SBATCH -n 64                                        # Number of cores
#SBATCH -p normal                                    # Partition to use
#SBATCH --mail-user=darmstrong4@islander.tamucc.edu  # Email for notifications
#SBATCH --mail-type=BEGIN,END,FAIL                   # Email on job start, end, and failure
#SBATCH --time=96:00:00                              # Maximum job runtime

# Define the path to the raw reads directory
RAW_READS="./00_fastq"  # Relative path from the repository root

# Define the output file for the checksum report
REPORT_FILE="./checksum_report.txt"

# Check if the raw reads directory exists
if [ ! -d "$RAW_READS" ]; then
  echo "Directory $RAW_READS does not exist. Exiting."
  exit 1
fi

# Clear or create the checksum report file
echo "Checksum Report" > "$REPORT_FILE"
echo "Generated on: $(date)" >> "$REPORT_FILE"
echo "------------------------------------------------------------" >> "$REPORT_FILE"

# Change to the raw reads directory
if cd "$RAW_READS"; then
  echo "Successfully changed to directory: $RAW_READS"
else
  echo "Failed to change directory to $RAW_READS. Exiting."
  exit 1
fi

# Verify the integrity of FASTQ files using md5 checksums
if ls *.md5 1> /dev/null 2>&1; then
  echo "Found checksum files. Starting verification..."
  for file in *.md5; do
    echo "Verifying checksum for: $file"
    {
      md5sum -c "$file" >> "../$REPORT_FILE" 2>&1
      if [ $? -eq 0 ]; then
        echo "Checksum verified successfully for $file." | tee -a "../$REPORT_FILE"
      else
        echo "Checksum verification failed for $file. Exiting." | tee -a "../$REPORT_FILE"
        exit 1
      fi
    }
  done
  echo "All checksums verified successfully!" | tee -a "../$REPORT_FILE"
else
  echo "No checksum files (.md5) found in $RAW_READS. Exiting." | tee -a "../$REPORT_FILE"
  exit 1
fi
