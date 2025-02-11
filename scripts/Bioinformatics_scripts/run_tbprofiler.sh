#!/bin/bash

# Directories to store output files
Reads="Fastq_files/"
SAMPLE_FILE="Sample_list.txt"

# Create necessary directories
mkdir -p bam vcf results

# To silence the citation option
parallel --citation --bibtex

# Run TB-Profiler in parallel
cat "$SAMPLE_FILE" | parallel --bar -j 2 tb-profiler profile -1 "$Reads"{}_1.fastq.gz -2 "$Reads"input/{}_2.fastq.gz -p {}
