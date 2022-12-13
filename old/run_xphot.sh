#!/bin/bash
#SBATCH -c 1                # Number of cores (-c)
#SBATCH -t 0-02:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p conroy_priority,test,shared,itc_cluster   # Partition to submit to
#SBATCH --mem=128000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o xmatch_phot.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e xmatch_phot.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --account=conroy_lab
module load python
source ~/.bashrc
conda activate outerhalo
cd /n/home03/vchandra/outerhalo/08_mage/
python -u 00_xmatch_phot.py
