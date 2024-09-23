#!/bin/bash

# Create necessary directories
mkdir -p Results/Trimmed Results/TrimmedUnpaired Results/Post_trimming_QC Results/sam Results/bam Results/vcf Results/bcf Results/sorted.bam Results/Aligned.stats Results/Multiqc_post_trimming

# Define a function to run Trimmomatic
trim_function() {
    infile=$1
    base=$(basename ${infile} _1.fastq.gz)
    Trimmed_R1=Results/Trimmed/${base}_1.trimmed.fastq.gz
    Trimmed_R2=Results/Trimmed/${base}_2.trimmed.fastq.gz
    TrimmedUn_R1=Results/TrimmedUnpaired/${base}_1un.trimmed.fastq.gz
    TrimmedUn_R2=Results/TrimmedUnpaired/${base}_2un.trimmed.fastq.gz
    trimmomatic PE -threads 20 -phred33 ${infile} reads/${base}_2.fastq.gz $Trimmed_R1 $TrimmedUn_R1 $Trimmed_R2 $TrimmedUn_R2 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True SLIDINGWINDOW:4:20 MINLEN:25
}

export -f trim_function

find reads -name '*_1.fastq.gz' | \
    parallel -j 3 trim_function {} {= s/_1.fastq.gz/_2.fastq.gz/ =}
