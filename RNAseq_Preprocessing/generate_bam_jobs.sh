#!/bin/bash

> bam_jobs.txt

mkdir -p bam_data


for sam in aligned_data/*.sam; do
    base=$(basename $sam .sam)

    echo "module purge; module load SAMtools; samtools view -@ 8 -Sb $sam | samtools sort -@ 8 -o bam_data/${base}_sorted.bam; samtools index bam_data/${base}_sorted.bam" >> bam_jobs.txt
done

