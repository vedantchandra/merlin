#!/bin/bash
#SBATCH -J r2023_04_30
#SBATCH -n 8 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 05:00:00 # Runtime
#SBATCH --mem-per-cpu 4500 # Memory
#SBATCH -p conroy_priority,shared,itc_cluster # Partition to submit to
#SBATCH --constraint='intel'
#SBATCH -o /n/holyscratch01/conroy_lab/vchandra/mage/logs/reduce/reduce_2023_04_30_v0.out
#SBATCH -e /n/holyscratch01/conroy_lab/vchandra/mage/logs/reduce/reduce_2023_04_30_v0.err

source activate pypeit2

cd /n/home03/vchandra/outerhalo/08_mage/pipeline/
echo 'CPU USED: ' 
cat /proc/cpuinfo | grep 'model name' | head -n 1
echo 'QUEUE NAME:' 
echo $SLURM_JOB_PARTITION
echo 'NODE NAME:' 
echo $SLURMD_NODENAME 

python -u radagast.py --dir=/n/holyscratch01/conroy_lab/vchandra/mage/data/2023_04_30/ --version=0  --skipred=False
