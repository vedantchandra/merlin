#!/bin/bash
#SBATCH -J rb2024_04_11
#SBATCH -n 8 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 05:00:00 # Runtime
#SBATCH --mem-per-cpu 4500 # Memory
#SBATCH -p conroy_priority,shared,itc_cluster # Partition to submit to
#SBATCH --constraint='intel'
#SBATCH -o /n/holystore01/LABS/conroy_lab/Lab/vchandra/mage/logs/reduce/reduce_b2024_04_11_v0.out
#SBATCH -e /n/holystore01/LABS/conroy_lab/Lab/vchandra/mage/logs/reduce/reduce_b2024_04_11_v0.err

source activate pypeit2

cd /n/home03/vchandra/outerhalo/08_mage/pipeline/
echo 'CPU USED: ' 
cat /proc/cpuinfo | grep 'model name' | head -n 1
echo 'QUEUE NAME:' 
echo $SLURM_JOB_PARTITION
echo 'NODE NAME:' 
echo $SLURMD_NODENAME 

python -u radagast.py --dir=/n/holystore01/LABS/conroy_lab/Lab/vchandra/mage/data/b2024_04_11/ --version=0  --skipred=False
