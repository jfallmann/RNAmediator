#!/bin/zsh
#SBATCH --job-name=RandRIssMed
#SBATCH --nodes=1
##SBATCH --nodelist=k72
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=96:00:00
#SBATCH --mem=20G
#SBATCH --mail-user=fall@bioinf.uni-leipzig.de
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --output=/scr/k70san2/fall/SLURM/RandRIssMed.%J.log
##SBATCH --error=/scr/k70san2/fall/SLURM/RandRIssMed.%J.err

echo "JOB = $SLURM_JOB_NAME"
echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"

ca rissmed
~/Projects/INPROGRESS/RIssmed/RIssmed/ConstraintPLFold.py -s GC$gc\.fa.gz -w 250 -l 150 -u $u -r raw -n unpaired -p paired -x sliding -o Random_$gc\_$u -z 10 -y $u --loglevel INFO -m 2
