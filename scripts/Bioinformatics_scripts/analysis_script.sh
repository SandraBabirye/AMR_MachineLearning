#!/bin/bash

# Set strict error handling
set -euo pipefail

# Check if the required arguments are provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <Sample_list.txt> <Results_output_directory>"
    exit 1
fi

# Assign arguments to variables
SAMPLE_FILE="$1"
RESULTS_DIR="$2"

# Define directories
Reads="${RESULTS_DIR}/Fastq_files/"
QC_Pre_trim="${RESULTS_DIR}/Pre_trim_QC"
MultiQC_pre_trim="${RESULTS_DIR}/Multiqc_Pre_trim"
Trimmed_reads="${RESULTS_DIR}/Trimmed"
Multiqc_post_trim="${RESULTS_DIR}/Multiqc_post_trimming"
QC_Post_trim="${RESULTS_DIR}/QC_Post_trim"
kraken_db_dir="Kraken2_database/minikraken2_v2_8GB_201904_UPDATE/"
kraken_reports_dir="${RESULTS_DIR}/Kraken"
bracken_reports_dir="${RESULTS_DIR}/Bracken"
classification_kraken_dir="${RESULTS_DIR}/classification_kraken"
SAM_DIR="${RESULTS_DIR}/sam"
BAM_DIR="${RESULTS_DIR}/bam"
SORTED_BAM_DIR="${RESULTS_DIR}/sorted_bam"
BCF_DIR="${RESULTS_DIR}/bcf"
VCF_DIR="${RESULTS_DIR}/vcf"
FILTERED_VCF_DIR="${RESULTS_DIR}/Filtered_vcfs"
SNPS_DIR="${RESULTS_DIR}/SNPS"
INDELS_DIR="${RESULTS_DIR}/INDELS"
POS_REF_ALT_DIR="${RESULTS_DIR}/POS_REF_ALT_Extracted"
REFERENCE_GENOME="ref/MTB_Ref.fasta"
THREADS=4

# Create necessary directories
mkdir -p "$Reads" "$QC_Pre_trim" "$MultiQC_pre_trim" "$Trimmed_reads" "$kraken_reports_dir" \
         "$bracken_reports_dir" "$classification_kraken_dir" "$SAM_DIR" "$BAM_DIR" "$SORTED_BAM_DIR" \
         "$BCF_DIR" "$VCF_DIR" "$FILTERED_VCF_DIR" "$SNPS_DIR" "$INDELS_DIR" "$POS_REF_ALT_DIR" "ref"

# Function to download samples
download_sample() {
    local sample=$1
    echo "Downloading $sample..."
    sratoolkit.3.0.7-ubuntu64/bin/fasterq-dump-orig.3.0.7 $sample --outdir "$Reads"
}
export -f download_sample

# Download samples in parallel
cat "$SAMPLE_FILE" | parallel -j $THREADS download_sample {}

# Run FastQC before trimming
find "$Reads" -name '*_1.fastq.gz' | parallel -j $THREADS fastqc -o "$QC_Pre_trim" --threads $THREADS {}

# Run MultiQC
multiqc -o "$MultiQC_pre_trim" "$QC_Pre_trim"

# Trimming function
trim_function() {
    infile="$1"
    base=$(basename "$infile" _1.fastq.gz)
    Trimmed_R1="$Trimmed_reads/${base}_1.trimmed.fastq.gz"
    Trimmed_R2="$Trimmed_reads/${base}_2.trimmed.fastq.gz"
    TrimmedUn_R1="$Trimmed_reads/${base}_1un.trimmed.fastq.gz"
    TrimmedUn_R2="$Trimmed_reads/${base}_2un.trimmed.fastq.gz"
    trimmomatic PE -threads $THREADS -phred33 "$infile" "${Reads}${base}_2.fastq.gz" \
        "$Trimmed_R1" "$TrimmedUn_R1" "$Trimmed_R2" "$TrimmedUn_R2" \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True SLIDINGWINDOW:4:20 MINLEN:25
}
export -f trim_function

# Run trimming in parallel
find "$Reads" -name '*_1.fastq.gz' | parallel -j $THREADS trim_function {}

# Run FastQC on trimmed reads
find "$Trimmed_reads" -name '*.fastq.gz' | parallel -j $THREADS fastqc -o "$QC_Post_trim" --threads $THREADS {}

# Run MultiQC after trimming
multiqc -o "$Multiqc_post_trim" "$QC_Post_trim"

# Kraken function
kraken_function() {
    trimmed_file_1=$1
    sample_name=$(basename "$trimmed_file_1" _1.trimmed.fastq.gz)
    trimmed_file_2="${Trimmed_reads}/${sample_name}_2.trimmed.fastq.gz"
    kraken2 --use-names --threads $THREADS --db "$kraken_db_dir" --fastq-input \
        --report "$kraken_reports_dir/${sample_name}.kraken_report" \
        --gzip-compressed --paired "$trimmed_file_1" "$trimmed_file_2" > "$classification_kraken_dir/${sample_name}.kraken"
    bracken -d "$kraken_db_dir" -i "$kraken_reports_dir/${sample_name}.kraken_report" \
        -l S -o "$bracken_reports_dir/${sample_name}.bracken"
}
export -f kraken_function

# Run Kraken2 classification in parallel
find "$Trimmed_reads" -name '*_1.trimmed.fastq.gz' | parallel -j 2 kraken_function {}

# Variant calling function
process_sample() {
    local sample=$1
    if [ ! -e "${REFERENCE_GENOME}.bwt" ]; then
        bwa index "$REFERENCE_GENOME"
    fi
    bwa mem -t $THREADS "$REFERENCE_GENOME" "${Trimmed_reads}/${sample}_1.trimmed.fastq.gz" \
        "${Trimmed_reads}/${sample}_2.trimmed.fastq.gz" > "$SAM_DIR/${sample}.sam"
        
    samtools view -S -b -o "$BAM_DIR/${sample}.bam" "$SAM_DIR/${sample}.sam"
    
    samtools sort "$BAM_DIR/${sample}.bam" -o "$SORTED_BAM_DIR/${sample}.sorted.bam"
    
    bcftools mpileup -O b -o "$BCF_DIR/${sample}.bcf" -f "$REFERENCE_GENOME" --threads $THREADS "$SORTED_BAM_DIR/${sample}.sorted.bam"
    
    bcftools call --ploidy 1 -m -v -o "$VCF_DIR/${sample}.vcf" "$BCF_DIR/${sample}.bcf"
    
    bcftools filter -i 'QUAL >= 30 && DP >= 10' "$VCF_DIR/${sample}.vcf" | \
        bcftools view -i 'FILTER="PASS"' -Oz -o "$FILTERED_VCF_DIR/${sample}_filtered.vcf.gz"
    
    bcftools view -v snps "$FILTERED_VCF_DIR/${sample}_filtered.vcf.gz" -Oz -o "$SNPS_DIR/${sample}_snps.vcf.gz"
    
    bcftools view -v indels "$FILTERED_VCF_DIR/${sample}_filtered.vcf.gz" -Oz -o "$INDELS_DIR/${sample}_indels.vcf.gz"
    
    bcftools query -f '%POS\t%REF\t%ALT\n' --print-header "$SNPS_DIR/${sample}_snps.vcf.gz" > "$POS_REF_ALT_DIR/${sample}.tsv"
    
    rm -f "$SAM_DIR/${sample}.sam" "$SORTED_BAM_DIR/${sample}.sorted.bam" "$BAM_DIR/${sample}.bam"
}
export -f process_sample

# Process samples in parallel
parallel -j $THREADS process_sample ::: $(cat "$SAMPLE_FILE")

echo "Pipeline completed successfully."
