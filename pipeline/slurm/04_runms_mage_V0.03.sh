#!/bin/bash
#SBATCH -J MS_MagE
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 05:00:00 # Runtime
#SBATCH --mem 5000 # Memory
#SBATCH -p conroy,shared,itc_cluster,serial_requeue # Partition to submit to
#SBATCH --constraint='intel'
#SBATCH -o /n/holyscratch01/conroy_lab/vchandra/mage/logs/mage/V0.03/%a.out
#SBATCH -e /n/holyscratch01/conroy_lab/vchandra/mage/logs/mage/V0.03/%a.err
#SBATCH --array=255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,275

source /n/home03/vchandra/warmup.sh

cd /n/home03/vchandra/outerhalo/08_mage/pipeline/
echo 'CPU USED: ' 
cat /proc/cpuinfo | grep 'model name' | head -n 1
echo 'QUEUE NAME:' 
echo $SLURM_JOB_PARTITION
echo 'NODE NAME:' 
echo $SLURMD_NODENAME 

python 04_runms_star.py --catalog=mage --ind=$SLURM_ARRAY_TASK_ID --version=V0.03 --npoints=500 --skipfit=0