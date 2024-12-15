# Soil Properties Drive Tissue-Specific Gene Expression in Miscanthus sinensis
Final project for CBB575 Group 8

Group members: Marlene Li, Ethan Burns, Pei-Wei Sun, and Zhenyang Zou (authors ordered by the workflow of the analyses)

The following serves as a documentation of the installation and configuration step. You can follow the instructions to reproduce the results.

## RNAseq_Preprocessing
This part contains the key files necessary for read cleaning, alignment, quantification, and normalization. Softwares used include Trimmomatic, HISAT2, SAMTool, featureCounts (from Subread), Python, and dSQ, all of which should be available via Yale HPC. `requirements.txt` for the conda environment is also provided. You can run the provided code files in the following order:

#### 1. Generate input file lists
By running:
```
bash generate_file_list.sh
```
You can get the `R1_files.txt` and `R2_files.txt`, which are a list of forward-read FASTQ files and a list of the reverse-read FASTQ files, respectively.

#### 2. Read cleaning via Trimmomatic
By running:
```
bash generate_trim_jobs.sh
```
You can get the `jobs.txt`, which contains jobs for trimming reads using Trimmomatic. After that, run:
```
bash run_dsq.sh
```
to submit tasks into SLURM job arrays for efficient parallelization.

#### 3. Mapping and alignment via HISAT2
By running:
```
bash generate_mapping_jobs.sh
```
You can get the `map-jobs.txt`, which contains mapping jobs. After that, run:
```
bash run_dsq2.sh
```
to submit tasks into SLURM job arrays for efficient parallelization.

#### 4. Convert to BAM format via SAMTool
By running:
```
bash generate_bam_jobs.sh
```
You can get the `bam_jobs.txt`, which contains jobs for converting .sam to .bam to prepare more efficient and practical files for featureCounts. After that, run:
```
bash run_dsq3.sh
```
to submit tasks into SLURM job arrays for efficient parallelization.

#### 5. Read quantification via featureCounts
By running:
```
bash generate_counts_jobs.sh
```
You can get the `fc_jobs.txt`, which contains jobs for counting features. After that, run:
```
bash run_dsq4.sh
```
to submit tasks into SLURM job arrays for efficient parallelization.

#### 6. Combine counts and calculate TPM
By running:
```
python generate_counts_and_tpm.py
```
You will get `counts_with_chr_and_length.csv` and `tpm_matrix.csv`, which will be used in the downstream analysis.

_Configuration Note:_

* Edit file paths and SLURM parameters (e.g., #SBATCH -p, #SBATCH -t) in .sh scripts to suit your environment.
* Edit the date in each `run_dsq*.sh` to the date you would like to execute the jobs.
* Ensure all tools are available in your PATH or loaded via module load.

## Differential Expression Analysis and Gene Ontology (GO) Analysis
Use additional modules R-bundle-CRAN/2024.06-foss-2022b on the Yale McCleary HPC R studio server to access packages. `DESeq2`, `ggrepel`, and `tidyverse` should all be available through there. R version R/4.4.1-foss-2022b was used. If needed, use the following code to install `DESeq2`.

```R
BiocManager::install("DESeq2")
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")

BiocManager::install("DESeq2")
```
## SNP Call and GWAS
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
