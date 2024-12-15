# Soil Properties Drive Tissue-Specific Gene Expression in Miscanthus sinensis
Final project for CBB575 Group 8
Group members: Marlene Li, Ethan Burns, Pei-Wei Sun, and Zhenyang Zou (authors ordered by the workflow of the analyses)

## SNP_call_and_GWAS
This part contains the key files necessary for identifying SNPs, conducting GWAS, mapping significant variants to genes, and performing enrichment analysis. Please follow the order below:

#### 1. `SNP Calling and Variant Generation.sh`
Run this first to call variants from BAM files using samtools and bcftools. It generates a filtered VCF suitable for GWAS.

#### 2. `GWAS Analysis with PLINK2.sh`
Next, run this script to perform GWAS for each trait using PLINK2, producing association statistics (p-values, effect sizes).

#### 3. `Identifying Significant Variants and Mapping to MisinGenes.sh`
Use this to extract significant variants from the GWAS results and map them to M. sinensis genes (MisinGenes) via bedtools.

#### 4. `new_gene.csv`
This file maps MisinGenes to their Arabidopsis thaliana orthologs, enabling functional annotation in the next step.

#### 5. `Converting MisinGenes and Running Enrichment.R`
Finally, run this R script to convert MisinGenes into Arabidopsis gene sets, perform GO/pathway enrichment with g:Profiler, and generate summary plots.
