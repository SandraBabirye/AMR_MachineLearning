#!/bin/bash

# Define paths and directories
fastqc_dir="FastQC_Results"
kraken_db_dir="minikraken2_v2_8GB_201904_UPDATE/"
trimmed_dir="Results/Trimmed/"
kraken_reports_dir="Kraken"
bracken_reports_dir="Bracken"
classification_kraken_dir="classification_kraken"
SAM_DIR="Results/sam"
BAM_DIR="Results/bam"
SORTED_BAM_DIR="Results/sorted.bam"
BCF_DIR="Results/bcf"
VCF_DIR="Results/vcf"
FILTERED_VCF_DIR="Results/Filtered_vcfs"
SNPS_DIR="Results/SNPS"
INDELS_DIR="Results/INDELS"
POS_REF_ALT_DIR="Results/POS_REF_ALT_Extracted"
PROJECT_READS_DIR="project_reads"
SAMPLE_FILE="Sample_list.txt"
REFERENCE_GENOME="MTB_Ref.fasta"
num_jobs=4

# Create necessary directories
mkdir -p "$fastqc_dir" "$kraken_reports_dir" "$bracken_reports_dir" "$classification_kraken_dir"
mkdir -p "$SAM_DIR" "$BAM_DIR" "$SORTED_BAM_DIR" "$BCF_DIR" "$VCF_DIR" "$FILTERED_VCF_DIR" "$SNPS_DIR" "$INDELS_DIR" "$POS_REF_ALT_DIR" "$PROJECT_READS_DIR"
mkdir -p Results/Trimmed Results/TrimmedUnpaired Results/Post_trimming_QC Results/Multiqc_post_trimming

# Step 1: Run FastQC on raw reads (before trimming)
echo "Running FastQC on raw reads..."
find reads -name '*_1.fastq.gz' | \
    parallel -j 4 fastqc -o "$fastqc_dir" --threads 4 {}

# Step 2: Perform trimming using Trimmomatic (in parallel)
trim_function() {
    infile=$1
    base=$(basename ${infile} _1.fastq.gz)
    Trimmed_R1="$trimmed_dir/${base}_1.trimmed.fastq.gz"
    Trimmed_R2="$trimmed_dir/${base}_2.trimmed.fastq.gz"
    TrimmedUn_R1="Results/TrimmedUnpaired/${base}_1un.trimmed.fastq.gz"
    TrimmedUn_R2="Results/TrimmedUnpaired/${base}_2un.trimmed.fastq.gz"
    trimmomatic PE -threads 20 -phred33 ${infile} reads/${base}_2.fastq.gz $Trimmed_R1 $TrimmedUn_R1 $Trimmed_R2 $TrimmedUn_R2 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True SLIDINGWINDOW:4:20 MINLEN:25

    ## Removing the unpaired reads
    rm -rf $TrimmedUn_R1
    rm -rf $TrimmedUn_R2

}

export -f trim_function

find reads -name '*_1.fastq.gz' | \
    parallel -j 3 trim_function {} {= s/_1.fastq.gz/_2.fastq.gz/ =}

# Step 3: Run FastQC on trimmed reads
echo "Running FastQC on trimmed reads..."
find "$trimmed_dir" -name '*_1.trimmed.fastq.gz' | \
    parallel -j 4 fastqc -o "$fastqc_dir" --threads 4 {}

# Step 4: Kraken and Bracken Taxonomy Investigation (in parallel)
echo "Running Kraken and Bracken analysis..."

kraken_function() {
    trimmed_file_1=$1
    sample_name=$(basename "$trimmed_file_1" _1.trimmed.fastq.gz)
    trimmed_file_2="${trimmed_dir}${sample_name}_2.trimmed.fastq.gz"
    
    kraken2 --use-names --threads 4 --db "$kraken_db_dir" --fastq-input --report "$kraken_reports_dir/${sample_name}.kraken" --gzip-compressed --paired "$trimmed_file_1" "$trimmed_file_2" > "${classification_kraken_dir}/${sample_name}.kraken"
    bracken -d "$kraken_db_dir" -i "$kraken_reports_dir/${sample_name}.kraken" -l S -o "$bracken_reports_dir/${sample_name}.bracken"
}

export -f kraken_function

find "$trimmed_dir" -name '*_1.trimmed.fastq.gz' | \
    parallel -j 4 kraken_function {}

# Step 5: Variant calling process (in parallel)
process_sample() {
    sample=$1
    
    # Indexing the reference genome (if not indexed already)
    if [ ! -e "$REFERENCE_GENOME.bwt" ]; then
        bwa index $REFERENCE_GENOME
    fi
    
    # Step 1: Alignment with BWA
    bwa mem -t 4 $REFERENCE_GENOME "$PROJECT_READS_DIR/${sample}_1.trimmed.fastq.gz" "$PROJECT_READS_DIR/${sample}_2.trimmed.fastq.gz" > "$SAM_DIR/${sample}.sam"

    # Step 2: Convert SAM to BAM and sort with Samtools
    samtools view -S -b -o "$BAM_DIR/${sample}.bam" "$SAM_DIR/${sample}.sam"
    samtools sort "$BAM_DIR/${sample}.bam" > "$SORTED_BAM_DIR/${sample}.sorted.bam"

    # Step 3: Variant calling with BCFtools
    bcftools mpileup -O b -o "$BCF_DIR/${sample}.bcf" -f $REFERENCE_GENOME --threads 8 "$SORTED_BAM_DIR/${sample}.sorted.bam"

    # Step 4: Convert BCF to VCF
    bcftools call --ploidy 1 -m -v -o "$VCF_DIR/${sample}.vcf" "$BCF_DIR/${sample}.bcf"

    # Step 5: Filtering the variants for Quality and Depth of coverage
    bcftools filter -i 'QUAL >= 30 && DP >= 10' "$VCF_DIR/${sample}.vcf" | bcftools view -i 'FILTER="PASS"' -Oz -o "$FILTERED_VCF_DIR/${sample}_filtered.vcf.gz"

    # Step 6: Extract SNPs from the VCF file
    bcftools view -v snps "$FILTERED_VCF_DIR/${sample}_filtered.vcf.gz" -Oz -o "$SNPS_DIR/${sample}_snps.vcf.gz"

    # Step 7: Extract Indels
    bcftools view -v indels "$FILTERED_VCF_DIR/${sample}_filtered.vcf.gz" -Oz -o "$INDELS_DIR/${sample}_indels.vcf.gz"

    # Step 8: Extract Position, Reference allele, and Alternative allele
    bcftools query -f '%POS\t%REF\t%ALT\n' --print-header "$SNPS_DIR/${sample}_snps.vcf.gz" > "$POS_REF_ALT_DIR/${sample}.tsv"

    # Step 9: Clean up intermediate files
    rm -rf "$SAM_DIR/${sample}.sam"
    rm -rf "$SORTED_BAM_DIR/${sample}.sorted.bam"
    rm -rf "$BAM_DIR/${sample}.bam"
}

export -f process_sample

# Run the variant calling in parallel for each input sample
mapfile -t input_files < "$SAMPLE_FILE"
parallel -j "$num_jobs" process_sample ::: "${input_files[@]}"

echo "Kraken, Bracken, and Variant calling analysis completed."

