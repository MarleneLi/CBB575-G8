#!/bin/bash

> jobs.txt


while read R1 R2; do
    base=$(basename $R1 _R1.fastq.gz)
    echo "module purge; module load Trimmomatic/0.39-Java-11; java -jar \$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -phred33 -threads 50 $R1 $R2 ${base}_R1_paired.fastq.gz ${base}_R1_unpaired.fastq.gz ${base}_R2_paired.fastq.gz ${base}_R2_unpaired.fastq.gz SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:36" >> jobs.txt
done < <(paste R1_files.txt R2_files.txt)

