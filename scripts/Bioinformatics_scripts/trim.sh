#!/bin/bash

# Define input reads directory
Reads="Fastq_files/"

# Create necessary directories
mkdir -p Results/Trimmed Results/TrimmedUnpaired

# Define a function to run Trimmomatic
trim_function() {
    infile="$1"
    base=$(basename "$infile" _1.fastq.gz)

    # Define output file paths
    Trimmed_R1="Results/Trimmed/${base}_1.trimmed.fastq.gz"
    Trimmed_R2="Results/Trimmed/${base}_2.trimmed.fastq.gz"
    TrimmedUn_R1="Results/TrimmedUnpaired/${base}_1un.trimmed.fastq.gz"
    TrimmedUn_R2="Results/TrimmedUnpaired/${base}_2un.trimmed.fastq.gz"

    # Run Trimmomatic
    trimmomatic PE -threads 20 -phred33 \
        "$infile" "${Reads}${base}_2.fastq.gz" \
        "$Trimmed_R1" "$TrimmedUn_R1" "$Trimmed_R2" "$TrimmedUn_R2" \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True SLIDINGWINDOW:4:20 MINLEN:25
}

export -f trim_function

# Run trimming in parallel
find "$Reads" -name '*_1.fastq.gz' | parallel -j 3 trim_function {}
