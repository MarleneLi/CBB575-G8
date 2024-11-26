module purge
module load dSQ

dsq --job-file bam_jobs.txt --mem-per-cpu 4g -t 20:00:00 --cpus-per-task 50 --mail-type ALL

sbatch dsq-bam_jobs-2024-11-21.sh
