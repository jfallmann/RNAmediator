#!/bin/zsh
#SBATCH --job-name=RissMEd
#SBATCH --nodes=1
#SBATCH --nodelist=k73
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=96:00:00
#SBATCH --mem=50G
#SBATCH --mail-user=fall@bioinf.uni-leipzig.de
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --output=/scr/k70san2/fall/SLURM/RIssMed.%J.log
##SBATCH --error=/scr/k70san2/fall/SLURM/RIssMed.%J.err

echo "JOB = $SLURM_JOB_NAME"
echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"

ca rissmed
~/Projects/INPROGRESS/RIssmed/RIssmed/ConstraintPLFold.py -s 65_Fu_targets.fa.gz -w 250 -l 150 -u 5 -g TargetGenes_65_Fu.bed -r raw -n unpaired -p paired -x 65_peaks_for_folding.bed -o 65_local_reproduce -z 30 --loglevel DEBUG -m 2
