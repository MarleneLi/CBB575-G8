#!/bin/bash


mkdir -p counts


> fc_jobs.txt

for bam in bam_data/*_sorted.bam; do
    base=$(basename $bam _sorted.bam)
    echo "module purge; module load Subread; featureCounts -T 40 -p -F GFF -t gene -g ID -a genome_data/Msinensis_497_v7.1.gene.gff3 -o counts/${base}_counts.txt $bam" >> fc_jobs.txt
done

