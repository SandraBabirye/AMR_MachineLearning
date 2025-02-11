#!/bin/bash

# Define paths and directories
kraken_db_dir="Kraken2_database/minikraken2_v2_8GB_201904_UPDATE/"
Trimmed_reads="Results/Trimmed/"
kraken_reports_dir="Results/Kraken"
bracken_reports_dir="Results/Bracken"
classification_kraken_dir="Results/classification_kraken"

# Create output directories if they don't exist
mkdir -p "$kraken_reports_dir" "$bracken_reports_dir" "$classification_kraken_dir"

# Define Kraken2 and Bracken classification function
kraken_function() {
    trimmed_file_1=$1
    sample_name=$(basename "$trimmed_file_1" _1.trimmed.fastq.gz)
    trimmed_file_2="${Trimmed_reads}${sample_name}_2.trimmed.fastq.gz"
    
    # Run Kraken2
    kraken2 --use-names --threads 4 --db "$kraken_db_dir" --fastq-input \
        --report "$kraken_reports_dir/${sample_name}.kraken_report" \
        --gzip-compressed --paired "$trimmed_file_1" "$trimmed_file_2" > "$classification_kraken_dir/${sample_name}.kraken"

    # Run Bracken (Species-level estimation)
    bracken -d "$kraken_db_dir" -i "$kraken_reports_dir/${sample_name}.kraken_report" \
        -l S -o "$bracken_reports_dir/${sample_name}.bracken"
}

export -f kraken_function

# Run classification in parallel
find "$Trimmed_reads" -name '*_1.trimmed.fastq.gz' | parallel -j 2 kraken_function {}
