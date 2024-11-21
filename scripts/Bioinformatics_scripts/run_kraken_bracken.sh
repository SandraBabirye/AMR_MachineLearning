#!/bin/bash

# Specify the Kraken database directory
kraken_db_dir="minikraken2_v2_8GB_201904_UPDATE/"

# Specify the directory containing the trimmed reads
trimmed_dir="Trimmed/"

# Create directories for Kraken and Bracken reports
mkdir -p Kraken Bracken classiffication_kraken

# Loop through each sample's paired trimmed read files in the directory
for trimmed_file_1 in "$trimmed_dir"*_1.trimmed.fastq.gz; do
    # Extract the sample name by removing the "_1.trimmed.fastq.gz" suffix
    sample_name=$(basename "$trimmed_file_1" _1.trimmed.fastq.gz)
    
    # Define the paths to the paired trimmed read files
    trimmed_file_2="${trimmed_dir}${sample_name}_2.trimmed.fastq.gz"
    
    # Run Kraken with the specified command
    kraken2 --use-names --threads 4 --db "$kraken_db_dir" --fastq-input --report "Kraken/${sample_name}.kraken"  --gzip-compressed --paired "$trimmed_file_1" "$trimmed_file_2" > "classiffication_kraken/{sample_name}.kraken" &
    
    # Wait for the Kraken process to finish before running Bracken
    wait $!
    
    # Run Bracken using the Kraken report
    bracken -d "$kraken_db_dir" -i "Kraken/${sample_name}.kraken" -l S -o "Bracken/${sample_name}.bracken" &
done

# Wait for all Bracken processes to finish
wait

echo "Kraken and Bracken analysis completed."
