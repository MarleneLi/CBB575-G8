#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --partition=day
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --job-name=snp_calling
#SBATCH -o /home/cbb575_zz469/Final/logs/snp_calling_%j.out
#SBATCH -e /home/cbb575_zz469/Final/logs/snp_calling_%j.err
# Directories and reference
BAM_DIR="/home/cbb575_zz469/Final/data/bam"
REF="/home/cbb575_zz469/Final/reference/Msinensis_497_v7.0.fa"
BAM_LIST="/home/cbb575_zz469/Final/data/bam/bam_list.txt"

# Output directory in scratch for large files
SCRATCH_DIR="/vast/palmer/scratch/cbb575/cbb575_zz469"

# Load modules (adjust versions as per your HPC environment)
module load GCCcore/10.2.0
module load SAMtools/1.16-GCCcore-10.2.0
module load HTSlib/1.16-GCCcore-10.2.0
module load BCFtools/1.16-GCCcore-10.2.0
module load VCFtools/0.1.16-GCCcore-10.2.0
module load parallel

# 1. Index each BAM file
for sp in $(cat $BAM_LIST); do
    samtools index -@ 4 ${BAM_DIR}/${sp}
done

#ls ${BAM_DIR}/*.bai


# 2. Create a chromosome list from the reference
samtools faidx $REF
cut -f1,2 ${REF}.fai > try
awk 'BEGIN{ FS=OFS="\t" } {$1 = $1 FS "0" }1' try > chromosomes.bed
rm try

# 3. Run bcftools mpileup/call for each chromosome in parallel
cat chromosomes.bed | awk '{print $1}' | parallel -j 2 "
bcftools mpileup -f $REF -r {} ${BAM_DIR}/*.bam | \
bcftools call -a GQ,GP -vm -o {}.vcf
"


# 4. Concatenate all chromosome VCF files into one merged VCF
bcftools concat *.vcf -O v -o ${SCRATCH_DIR}/merged.vcf

# Count the number of SNPs
bcftools view -v snps ${SCRATCH_DIR}/merged.vcf | wc -l

# 5. Filter SNPs using vcftools and store result in scratch
vcftools --vcf ${SCRATCH_DIR}/merged.vcf --max-missing 1 --recode --recode-INFO-all --out ${SCRATCH_DIR}/filt

# Final filtered file is ${SCRATCH_DIR}/filt.recode.vcf

