#!/bin/bash
#SBATCH -J fitMAGE
#SBATCH -c 1                # Number of cores (-c)
#SBATCH -t 0-03:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p conroy_priority,itc_cluster,shared,serial_requeue   # Partition to submit to
#SBATCH --constraint="intel"
#SBATCH --mem-per-cpu=3500           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o logs/msfit_%a.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e logs/msfit_%a.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --account=conroy_lab
#SBATCH --array=0-79

module load python
source /n/home03/vchandra/.bashrc
source activate outerhalo
cd /n/home03/vchandra/outerhalo/08_mage/

python -u 01_runstar.py "${SLURM_ARRAY_TASK_ID}" --version='h3'
