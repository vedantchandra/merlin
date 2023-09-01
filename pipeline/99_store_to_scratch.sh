#! /bin/sh
#SBATCH -c 1                # Number of cores (-c)
#SBATCH -t 0-06:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p conroy_priority,shared,itc_cluster   # Partition to submit to
#SBATCH --mem=16000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o slurm/backup.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e slurm/backup.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --account=conroy_lab

cp -R /n/holystore01/LABS/conroy_lab/Lab/vchandra/mage /n/holyscratch01/conroy_lab/vchandra/
