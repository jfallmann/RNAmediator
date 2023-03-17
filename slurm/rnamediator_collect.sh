#!/bin/bash
#SBATCH --job-name=RNAmedRan
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=196:00:00
#SBATCH --mem=20G
#SBATCH --mail-user=fall@bioinf.uni-leipzig.de
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --output=RNAmedRan.%J.log
##SBATCH --error=RNAmedRan.%J.err

echo "JOB = $SLURM_JOB_NAME"
echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"

RNAmediator_collect_plfold -u 9 -g ${GROUP}.bed.gz -b 1 -c 0 -p ${WINDOW},${SPAN} -o Collection_${GROUP} -d ${GROUP}_${WINDOW}_${SPAN}_fold