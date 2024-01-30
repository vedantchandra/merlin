#!/bin/bash
#SBATCH -J MS_MagE
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 15:00:00 # Runtime
#SBATCH --mem 5000 # Memory
#SBATCH -p conroy,shared,itc_cluster,serial_requeue # Partition to submit to
#SBATCH --constraint='intel'
#SBATCH -o /n/holyscratch01/conroy_lab/vchandra/mage/logs/mage/V0.05/%a.out
#SBATCH -e /n/holyscratch01/conroy_lab/vchandra/mage/logs/mage/V0.05/%a.err
#SBATCH --array=63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,229,270,271,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,289,290,291,292,293,294,295,296,297,298,299,300,301

source activate outerhalo

cd /n/home03/vchandra/outerhalo/08_mage/pipeline/
echo 'CPU USED: ' 
cat /proc/cpuinfo | grep 'model name' | head -n 1
echo 'QUEUE NAME:' 
echo $SLURM_JOB_PARTITION
echo 'NODE NAME:' 
echo $SLURMD_NODENAME 

python 04_runms_star.py --catalog=mage --ind=$SLURM_ARRAY_TASK_ID --version=V0.05 --npoints=500 --skipfit=0
