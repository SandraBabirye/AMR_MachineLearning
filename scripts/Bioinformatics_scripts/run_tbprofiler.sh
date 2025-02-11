#!/bin/bash

#Directories to store output files

mkdir bam vcf results

#To silence the citation option

parallel --citation --bibtex

cat Samples_list.txt | parallel --bar -j 2 tb-profiler profile -1 $Reads{}_1.fastq.gz -2  Research_container/input/{}_2.fastq.gz -p {} 
