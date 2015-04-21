#SBATCH -J blast # A single job name for the array

#SBATCH -n 2 # Number of cores

#SBATCH -N 1 # All cores on one machine

#SBATCH -p serial_requeue # Partition

#SBATCH --mem 4000 # Memory request

#SBATCH -t 0-2:00 # 2 hours (D-HH:MM)

#SBATCH -o blast_%A_%a.out # Standard output

#SBATCH -e blast_%A_%a.err # Standard error



blastp -num_threads 2-query input${SLURM_ARRAY_TASK_ID}.fasta -db refactor_test -out output${SLURM_ARRAY_TASK_ID}.blast