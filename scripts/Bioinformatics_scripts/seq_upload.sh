#!/bin/bash

# Create the project_reads directory if it doesn't exist
mkdir -p project_reads

# Path to the text file containing sample names
SAMPLE_FILE="samples.txt"

# Function to download a sample using fasterq-dump
download_sample() {
    local sample=$1
    echo "Downloading $sample..."
    sratoolkit.3.0.7-ubuntu64/bin/fasterq-dump-orig.3.0.7 $sample --outdir project_reads/
}

export -f download_sample

# Read the sample names from the file and run downloads in parallel
cat $SAMPLE_FILE | parallel -j 4 download_sample
