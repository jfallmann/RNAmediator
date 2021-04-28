#!/bin/zsh
#SBATCH --job-name=RIssMed
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
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
cd /scr/k70san2/fall/Constraints/miRNAs/mirwalk_coop
~/Projects/INPROGRESS/RIssmed/RIssmed/CollectConsResults.py -u 7 -z 15 -g ../targetgenes.bed.gz -b-1,1 -c 0 -p 250,150 -o Collection_mirnas
