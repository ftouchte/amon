#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --array=0-834
#SBATCH --job-name=ftc_kf_cooking_%a
#SBATCH --output=/farm_out/%u/%x-%A-%a.out
#SBATCH --error=/farm_out/%u/%x-%A-%a.err
#SBATCH --time=07:00:00
#SBATCH --ntasks=1
#SBATCH --mem=2000

run=22712
i=$SLURM_ARRAY_TASK_ID

echo "Running job $i"

/w/hallb-scshelf2102/clas12/users/touchte/data/simu/run-job.sh $run $i
