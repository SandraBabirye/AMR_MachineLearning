## Machine Learning-based prediction of antibiotic resistance in Mycobacterium tuberculosis clinical isolates from Uganda 
Antimircobial resistance (AMR) is a growing challenge recognised as global health emergency posing a great threat to modern medicine. 
In this study we used both bioinformatics and machine learning approaches for analysis.

## Bioinformatics analysis for the whole genome sequence data
This analysis was performed using various key steps and tools/softwares in the command line.

## Overview of the steps
- Quality control
   1.Quality assessment using FastQC and MultiQC
   2.Trimming the reads using Trimmomatic
   3.Read taxonomy investigation
- Indexing the reference genome
- Read mapping
- Variant calling
- Variant filtering
- Variant annotation

## Pre-Requisites

Before running the script, ensure the following **tools** are installed together with their **dependencies**

- **FastQC**
- **Trimmomatic**
- **Kraken2 Database**: Download the Kraken2 database, such as `minikraken2_v2_8GB_201904_UPDATE`.
- **BWA**: for reference genome alignment.
- **Samtools**
- **BCFtools**





###### Bioinformatics analysis workflow
![Alt text](Figures/Bioinformatics_analysis_workflow.jpg)
#### Genomic data process

### Machine learning analysis workflow
In the current study, we ten machine learning algorithms including  **Multilayer Perceptron (MLP), Extreme Gradient Boosting classifier (XGBoost), Gradient Boosting (GBC), Adaptive Boosting (adaboost), logistic regression (LR), support vector machine (SVM), random forest (RF), Decision Trees (DT), Extra Trees Classifier (ETC), CatBoost** for the prediction of AMR for the antibiotics rifampicin (RIF), isoniazid (INH), ethambutol (EMB), and streptomycin (STM) based on WGS and clinical data (Age, sex, HIV status)  with one hot encoding.

This repository contains the required python scripts and associated data to train and test Machine Learning (ML) models using most classifiers supported by the python's ski-learn, and associated Bioinformatics packages in a Linux environment. The sample ML script for only one set of drug(rifampicin) has been provided as jupyter notebook file (.ipynb) to enhance users ability to break the scripts down into manageable and understandable sections.
![Alt text](Figures/ML_analysis_workflow.png)

### Please cite our article here
Babirye, S.R., Nsubuga, M., Mboowa, G. et al. Machine learning-based prediction of antibiotic resistance in Mycobacterium tuberculosis clinical isolates from Uganda. BMC Infect Dis 24, 1391 (2024). https://doi.org/10.1186/s12879-024-10282-7
