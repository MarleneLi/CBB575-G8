#!/bin/bash

cd /home/cbb575_cl2658/palmer_scratch/final_project/processed_data


> map-jobs.txt

for R1 in trim_data/*_R1_paired.fastq.gz; do
    R2=${R1/_R1_/_R2_}  
    base=$(basename $R1 _R1_paired.fastq.gz)

    echo "module purge; module load HISAT2; hisat2 -p 32 -x genome_data/Msinensis_index -1 $R1 -2 $R2 -S aligned_data/${base}.sam" >> map-jobs.txt
done

