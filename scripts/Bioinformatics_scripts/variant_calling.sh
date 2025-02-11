#!/bin/bash

# Define directories
Trimmed_reads="Results/Trimmed/"
SAM_DIR="Results/sam"
BAM_DIR="Results/bam"
SORTED_BAM_DIR="Results/sorted_bam"
BCF_DIR="Results/bcf"
VCF_DIR="Results/vcf"
FILTERED_VCF_DIR="Results/Filtered_vcfs"
SNPS_DIR="Results/SNPS"
INDELS_DIR="Results/INDELS"
POS_REF_ALT_DIR="Results/POS_REF_ALT_Extracted"
REFERENCE_GENOME="ref/MTB_Ref.fasta"

# Create necessary directories
mkdir -p "$SAM_DIR" "$BAM_DIR" "$SORTED_BAM_DIR" "$BCF_DIR" "$VCF_DIR" \
         "$FILTERED_VCF_DIR" "$SNPS_DIR" "$INDELS_DIR" "$POS_REF_ALT_DIR" "ref"

# Set the number of parallel jobs
num_jobs=4

# Read input sample names from the text file into an array
mapfile -t input_files < Sample_list.txt

# Define the function to process each sample in parallel
process_sample() {
    local sample=$1

    # Indexing the reference genome (only if not already indexed)
    if [ ! -e "${REFERENCE_GENOME}.bwt" ]; then
        bwa index "$REFERENCE_GENOME"
    fi

    # Step 1: Alignment with BWA
    bwa mem -t 4 "$REFERENCE_GENOME" \
        "${Trimmed_reads}${sample}_1.trimmed.fastq.gz" \
        "${Trimmed_reads}${sample}_2.trimmed.fastq.gz" > "${SAM_DIR}/${sample}.sam"

    # Step 2: Convert SAM to BAM
    samtools view -S -b -o "${BAM_DIR}/${sample}.bam" "${SAM_DIR}/${sample}.sam"

    # Step 3: Sort BAM file
    samtools sort "${BAM_DIR}/${sample}.bam" -o "${SORTED_BAM_DIR}/${sample}.sorted.bam"

    # Step 4: Variant calling with BCFtools
    bcftools mpileup -O b -o "${BCF_DIR}/${sample}.bcf" -f "$REFERENCE_GENOME" --threads 8 "${SORTED_BAM_DIR}/${sample}.sorted.bam"

    # Step 5: Convert BCF to VCF
    bcftools call --ploidy 1 -m -v -o "${VCF_DIR}/${sample}.vcf" "${BCF_DIR}/${sample}.bcf"

    # Step 6: Filter variants (QUAL ≥ 30 & DP ≥ 10)
    bcftools filter -i 'QUAL >= 30 && DP >= 10' "${VCF_DIR}/${sample}.vcf" | \
        bcftools view -i 'FILTER="PASS"' -Oz -o "${FILTERED_VCF_DIR}/${sample}_filtered.vcf.gz"

    # Step 7: Extract SNPs
    bcftools view -v snps "${FILTERED_VCF_DIR}/${sample}_filtered.vcf.gz" -Oz -o "${SNPS_DIR}/${sample}_snps.vcf.gz"

    # Step 8: Extract Indels
    bcftools view -v indels "${FILTERED_VCF_DIR}/${sample}_filtered.vcf.gz" -Oz -o "${INDELS_DIR}/${sample}_indels.vcf.gz"

    # Step 9: Extract Position, Reference Allele, and Alternative Allele
    bcftools query -f '%POS\t%REF\t%ALT\n' --print-header "${SNPS_DIR}/${sample}_snps.vcf.gz" > "${POS_REF_ALT_DIR}/${sample}.tsv"

    # Clean up intermediate files to save space
    rm -f "${SAM_DIR}/${sample}.sam" "${SORTED_BAM_DIR}/${sample}.sorted.bam" "${BAM_DIR}/${sample}.bam"
}

export -f process_sample

# Run the variant calling in parallel for each input sample
parallel -j "$num_jobs" process_sample ::: "${input_files[@]}"
