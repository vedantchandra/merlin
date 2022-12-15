#!/bin/bash
#SBATCH -c 1                # Number of cores (-c)
#SBATCH -t 0-01:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p conroy_priority,shared,itc_cluster   # Partition to submit to
#SBATCH --mem=128000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o /n/holyscratch01/conroy_lab/vchandra/mage/logs/xgall/xgall_%a.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e //n/holyscratch01/conroy_lab/vchandra/mage/logs/xgall/xgall_%a.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --account=conroy_lab
#SBATCH --array=0-26

module load python
source ~/.bashrc
conda activate outerhalo

cd /n/home03/vchandra/outerhalo/08_mage/pipeline/
python -u 02_xmatch_gall.py "${SLURM_ARRAY_TASK_ID}"