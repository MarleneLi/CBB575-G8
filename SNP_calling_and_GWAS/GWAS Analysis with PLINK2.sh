#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --partition=day
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --job-name=gwas_all_metals
#SBATCH -o /home/cbb575_zz469/Final/logs/gwas_all_metals_%j.out
#SBATCH -e /home/cbb575_zz469/Final/logs/gwas_all_metals_%j.err

module load PLINK/2_avx2_20221024

GWAS_DIR="/home/cbb575_zz469/Final"
SCRATCH_DIR="/vast/palmer/scratch/cbb575/cbb575_zz469"

traits=("Al" "Ca" "Fe" "K" "Cr" "Pb" "Mn" "As" "Cu" "Ni" "Zn" "S")

for t in "${traits[@]}"; do
    PHENO_FILE="phenotype_${t}.txt"
    if [ -f "$PHENO_FILE" ]; then
        plink2 --bfile ${GWAS_DIR}/gwas_data_biallelic \
               --pheno $PHENO_FILE \
               --pheno-col-nums 3 \
               --glm allow-no-covars \
               --out ${SCRATCH_DIR}/gwas_${t}
    else
        echo "Phenotype file for $t not found: $PHENO_FILE"
    fi
done
