## Machine Learning-based prediction of antibiotic resistance in Mycobacterium tuberculosis clinical isolates from Uganda 
Antimicrobial resistance (AMR) is a growing challenge recognised as global health emergency posing a great threat to modern medicine. 

This repository contains bioinformatics scripts, Jupyter notebooks, and data for training and testing machine learning (ML) algorithms implemented using the Scikit-learn library in Python. A sample ML script for predicting drug resistance to Rifampicin has been provided as a Jupyter notebook file (.ipynb). This format allows users to break down the script into manageable and understandable sections.

## Bioinformatics analysis for the whole genome sequence data
This analysis was performed using various key steps and tools/softwares in the command line.

#### Bioinformatics analysis workflow
![Alt text](Figures/Bioinformatics_analysis_workflow.jpg)

## Overview of the steps
- Quality control
  1.Quality assessment 
  2.Trimming the reads 
  3.Read taxonomy investigation
- Indexing the reference genome
- Read mapping
- Variant calling
- Variant filtering
- Variant annotation
- Lineage and drug resistance prediction 

**Note** We used the `snippy pipeline` that does variant calling and also generates a core genome alignment file for phylogenetic tree construction. However in this document we have provided detailed step by step information for some of the individual tools used in the pipeline. We have also created a bash script with all these steps merged into one.

## Pre-Requisites

Before running the script, ensure the following **tools** are installed together with their **dependencies**. These tools can be installed using Conda

- `FastQC`
- `Fasterqdump` from the SRAtool kit
- `Trimmomatic`
- `BWA`
- `Samtools`
- `BCFtools`
- `TBProfiler`
- `Snippy`
- `SnpEff`
- `Kraken2`
- `bracken`
- `parallel`

## Environment Setup

To ensure that you have all the necessary dependencies and tools installed , we have provided an environment YAML file in the `Envs` directory. Follow these steps to set up the environment:

1. Navigate to the `Envs` directory where the YAML file is located:
   ```
   cd Envs
   ```
   
2. Create and activate the environment using Conda (make sure Conda is installed on your system):
   ```
   conda env create -f environment.yml
   
   ```
4. Activate the environment. Replace <environment_name> with the actual environment name specified in the environment.yml file.
   ```
   conda activate <environment_name>

   ```
   
### File structure

- `Fastq_files`: Directory containing raw sequencing data files `(.fastq.gz)`.
- `Results/` :  Directory for all output results including trimmed reads, QC reports, and analysis results.
- `Kraken2 database` : Kraken database directory such as `minikraken2_v2_8GB_201904_UPDATE`
- `ref` : Directory containing the **reference genome**
- `Sample_list.txt` : A file containing the list of accession numbers for various samples to be dowloaded from a public repository for analysis.


### Detailed Pipeline Steps 

#### Setp 1: Download the raw reads from the public database

**Input**: `Sample_list.txt`

```
Reads="Fastq_files/"

# Create the project_reads directory if it doesn't exist
mkdir -p "$Reads"

# Path to the text file containing sample names
SAMPLE_FILE="Sample_list.txt"

# Function to download a sample using fasterq-dump
download_sample() {
    local sample=$1
    echo "Downloading $sample..."
    sratoolkit.3.0.7-ubuntu64/bin/fasterq-dump-orig.3.0.7 $sample --outdir "$Reads"
}

export -f download_sample

# Read the sample names from the file and run downloads in parallel
cat $SAMPLE_FILE | parallel -j 4 download_sample

```
#### Step 2: Run FastQC on Raw Reads (Before Trimming)

FastQC helps to assess the quality of the raw reads. MultiQC helps aggregate the single qc reposrts generated by FastQC to form one general qc report for all samples.

**Input**: Raw paired-end fastq files (_1.fastq.gz).

**Output**: Quality control reports saved in the FastQC_Results directory.

```
# Define the directory variable
QC_pre_trim="Results/Pre_trim_QC" 
MultiQC_pre_trim="Results/Multiqc_Pre_trim"

# Creating directories 
mkdir -p  "$QC_Pre_trim" "$MultiQC_pre_trim"

# Run FastQC in parallel on all *_1.fastq.gz files in the 'reads' directory
find reads -name '*_1.fastq.gz' | parallel -j 4 fastqc -o "$QC_Pre_trim" --threads 4 {}

# Run MultiQC
multiqc -o "$MultiQC_pre_trim" "$QC_Pre_trim"

```

### Step 3: Perform Trimming with Trimmomatic

Trimmomatic is used to trim low-quality sequences and adapter contamination from reads. The process is done in parallel to speed up the workflow.
**NB**: The arguments for trimmomatic can varry depending on the quality of your reads ie ```LEADING, TRAILING, MINLEN, SLIDINGWINDOW```etc. Also we used `TruSeq3-PE.fa` as the adapater sequence in this study. It also varry.

**Input**: Paired-end raw fastq files.

**Output**: Trimmed paired-end reads saved in the `Results/Trimmed/` directory and unpaired reads in the `Results/TrimmedUnpaired/` directory.
```
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
```

### Step 4: Run FastQC on Trimmed reads

After trimming, FastQC is run again to assess the quality of the trimmed reads.

**Input**: Trimmed paired-end fastq files.
**Output**: Quality control reports for trimmed reads.

```

# Define the directory variable
Trimmed_reads= "Results/Trimmed"
Multiqc_post_trim= "Results/Multiqc_post_trimming"
QC_Post_trim="Results/QC_Post_trim"

# Creating directories 
mkdir -p "$fastqc_dir" "$Multiqc_dir"

# Run FastQC on trimmed reads
find "$Trimmed_reads" -name '*.fastq.gz' | parallel -j 4 fastqc -o "$fastqc_dir" --threads 4 {}

# Run MultiQC to aggregate FastQC results
multiqc -o "$Multiqc_post_trim" "$QC_Post_trim"

```

### Step 5: Kraken and Bracken Taxonomy Classification
Kraken2 is used for taxonomic classification of sequencing reads. Bracken is then used to improve the abundance estimation of taxa at different ranks (e.g., species, genus).

**Input**: Trimmed paired-end fastq files.

**Output**: Kraken classification reports and Bracken abundance estimates.

```
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

```

### Step 6: Variant calling, Variant filtering, Variant annotation and VCF File manipulation

The variant calling process is performed in multiple steps for each sample:

- Alignment: Reads are aligned to the reference genome using BWA.
- SAM to BAM Conversion: SAM files are converted to BAM format using Samtools.
- Sorting: BAM files are sorted.
- Variant Calling: Variants are called using BCFtools.
- VCF Filtering: Variants are filtered based on quality and depth of coverage.
- Extracting SNPs and Indels: SNPs and Indels are extracted from the filtered VCF files.
- Extracting Positions, Reference, and Alternative Alleles: Position and alleles are extracted for downstream analysis.


**Input**: Trimmed paired-end fastq files.

**Output**: Various output files including aligned BAM files, sorted BAM files, VCF files, filtered VCF files, and extracted SNP and Indel files.

```
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

```
**Note** : To improve efficiency, most tasks (FastQC, trimming, Kraken classification, and variant calling) are run in parallel using GNU Parallel. This allows the pipeline to process multiple samples simultaneously, speeding up the analysis for large datasets.

### Step 7: Running the bash script
When running the final script provide the path to the sample list as ` `Sample_list.txt`  and that for the output results generated during the analysis `Results`

```
bash analysis_script.sh  Sample_list.txt Results
```

### Additional step specific to MTB isolates which i used in this project: Lineage and drug resistance prediction using TBProfiler
```
# Directories to store output files
Reads="Fastq_files/"
SAMPLE_FILE="Sample_list.txt"

# Create necessary directories
mkdir -p bam vcf results

# To silence the citation option
parallel --citation --bibtex

# Run TB-Profiler in parallel
cat "$SAMPLE_FILE" | parallel --bar -j 2 tb-profiler profile -1 "$Reads"{}_1.fastq.gz -2 "$Reads"input/{}_2.fastq.gz -p {}

```


### Machine learning (ML) analysis

In the current study, we trained ten machine learning algorithms to predict antimicrobial resistance (AMR) for the antibiotics Rifampicin (RIF), Isoniazid (INH), Ethambutol (EMB), and Streptomycin (STM). The algorithms used include:

1. **Extreme Gradient Boosting classifier (XGBoost)**
2. **Gradient Boosting (GBC)**
3. **Adaptive Boosting (adaboost)**
4. **Logistic Regression (LR)**
5. **Support vector machine (SVM)**
6. **Random Forest (RF)**
7. **Decision Trees (DT)**
8. **Extra Trees Classifier (ETC)**
9. **CatBoost**
10. **Multilayer Perceptron (MLP)**
    
The prediction was based on whole-genome sequencing (WGS) and clinical data (age, sex, and HIV status), with categorical features encoded using one-hot encoding.

**Note**: The SNPmatrix which consisted of genomic variations across the entire genome and the clinica data variables were the  the predictor variables and the phenotypic drug susceptibility testing (DST) data for the drugs was the target variable.


#### ML analysis workflow
![Alt text](Figures/ML_analysis_workflow.png)

### Please cite our article here
Babirye, S.R., Nsubuga, M., Mboowa, G. et al. Machine learning-based prediction of antibiotic resistance in Mycobacterium tuberculosis clinical isolates from Uganda. BMC Infect Dis 24, 1391 (2024). https://doi.org/10.1186/s12879-024-10282-7
