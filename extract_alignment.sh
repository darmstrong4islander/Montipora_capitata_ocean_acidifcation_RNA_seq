#!/bin/bash

# Script to extract alignment percentages from the log file

# Define the error file
ERROR_FILE="mcap_rna_cnt.err"

# Define the output file
OUTPUT_FILE="alignment_percentages.txt"

# Extract alignment rates and corresponding sample names
grep -E "overall alignment rate|fastp --in1" "$ERROR_FILE" | \
awk '
    /fastp --in1/ { 
        split($3, arr, "/"); 
        sample=arr[3]; 
        sub("_R1_001.fastq.gz", "", sample)
    }
    /overall alignment rate/ { 
        print sample, $1, $2, $3, $4, $5
    }
' > "$OUTPUT_FILE"

echo "Alignment percentages saved to $OUTPUT_FILE"

