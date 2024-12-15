#!/bin/bash
#SBATCH --job-name=extract_misingenes
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -o extract_misingenes_%j.out
#SBATCH -e extract_misingenes_%j.err

module load parallel
# Load any other modules you need, if required

cd /vast/palmer/scratch/cbb575/cbb575_zz469  # Adjust to your working directory

traits=("Al" "Ca" "Fe" "K" "Cr" "Pb" "Mn" "As" "Cu" "Ni" "Zn" "S")
p_threshold=1e-4

for t in "${traits[@]}"; do
    gwas_file="gwas_${t}.${t}.glm.linear"
    # Extract significant variants
    awk -v p=$p_threshold 'NR>1 && $12 < p {print $1":"$2}' $gwas_file > ${t}_hits.txt
    
    # Grep these variants in variants_with_genes_id.bed
    grep -F -f ${t}_hits.txt variants_with_genes_id.bed | awk '{print $8}' | sort | uniq > ${t}_misingenes.txt
done
