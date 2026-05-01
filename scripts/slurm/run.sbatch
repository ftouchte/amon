#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --job-name=rgl_ahdc_wire_position
#SBATCH --output=/farm_out/%u/%x-%A-%a.out
#SBATCH --error=/farm_out/%u/%x-%A-%a.err
#SBATCH --time=07:00:00
#SBATCH --ntasks=1
#SBATCH --mem=6000

# Files location
workDir="/lustre24/expphy/volatile/clas12/touchte/alignment-cooking"
input_list="$workDir/input_list.txt"
output_list="$workDir/output_list.txt"

# Read file stored in line number slurm_id + 1 the relevant files
decfile=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p"  "$input_list")
recfile=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p"  "$output_list")

export CCDB_CONNECTION="sqlite:////work/clas12/users/touchte/ccdb_2026-04-26.sqlite"

coatjava="/w/hallb-scshelf2102/clas12/users/touchte/coatjava"

yamlfile="/lustre24/expphy/volatile/clas12/touchte/alignment-cooking/config.yaml"


# Summary
echo "CCDB_CONNECTION : $CCDB_CONNECTION"
echo "coatjava : $coatjava"
echo "yamlfile : $yamlfile"
echo "command  : $coatjava/coatjava/bin/recon-util -i $decfile  -o $recfile -y $yamlfile"

# Check CCDB
echo "> CCDB table check"
ccdb -r 22712 dump /calibration/alert/ahdc/time_offsets | grep "192 *0\.00"

# Uncomment to run command
#$coatjava/coatjava/bin/recon-util -i "$decfile"  -o "$recfile "-y "$yamlfile"
