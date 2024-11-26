module purge
module load dSQ

dsq --job-file jobs.txt --mem-per-cpu 4g -t 20:00:00 --cpus-per-task 50 --mail-type ALL

sbatch dsq-jobs-2024-11-19.sh
