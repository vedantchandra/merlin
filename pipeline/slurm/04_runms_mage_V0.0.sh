#!/bin/bash
#SBATCH -J MS_MagE
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 05:00:00 # Runtime
#SBATCH --mem 5000 # Memory
#SBATCH -p conroy,shared,itc_cluster,serial_requeue # Partition to submit to
#SBATCH --constraint='intel'
#SBATCH -o /n/holyscratch01/conroy_lab/vchandra/mage/logs/mage/V0.0/%a.out
#SBATCH -e /n/holyscratch01/conroy_lab/vchandra/mage/logs/mage/V0.0/%a.err
#SBATCH --array=0,1,3,4,6,7,8,9,10,11,12,13,14,15,16,17,21,23,25,26,27,28,29,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,53,54,55,56,57,60,61,62,63,64,65,66,72,73,75,78,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,126,127,128,129,130,131,132,133,134,135,136,137,138,139,142,143,146,147,148,149,150,151,152,153,154,155,156,166,168,172,173,175,176,182,183,184,185,186,187,188,190,191,193,195,196,197,199,244,257,258,268,269,270,271,272,273,274,275

source /n/home03/vchandra/warmup.sh

cd /n/home03/vchandra/outerhalo/08_mage/pipeline/
echo 'CPU USED: ' 
cat /proc/cpuinfo | grep 'model name' | head -n 1
echo 'QUEUE NAME:' 
echo $SLURM_JOB_PARTITION
echo 'NODE NAME:' 
echo $SLURMD_NODENAME 

python 04_runms_star.py --catalog=mage --ind=$SLURM_ARRAY_TASK_ID --version=V0.0 --npoints=500 --skipfit=0