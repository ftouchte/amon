#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --job-name=rgl_ahdc_new_translation_table
#SBATCH --output=/farm_out/%u/%x-%A-%a.out
#SBATCH --error=/farm_out/%u/%x-%A-%a.err
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --mem=6000

# Files location
workDir="/lustre24/expphy/volatile/clas12/touchte/new-translation-table" # should match the one in submit_slurm_job.sbatch

# The following files are created in submit_slurm_job.sbatch
evio_list="$workDir/evio_list.txt"
dec_list="$workDir/dec_list.txt"
rec_list="$workDir/rec_list.txt"

# Read file stored in line number slurm_id + 1 the relevant files
eviofile=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p"  "$evio_list")
decfile=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p"  "$dec_list")
recfile=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p"  "$rec_list")

# Place to store temporary files
export TMPDIR="/scratch/touchte/tmpfs"

# Set path to the CCDB
export CCDB_CONNECTION="sqlite:////work/clas12/users/touchte/ccdb_2026-04-26.sqlite"

# path to coatjava
coatjava="/w/hallb-scshelf2102/clas12/users/touchte/coatjava"

# Path to the YAML file
yamlfile="/volatile/clas12/touchte/new-translation-table/config.yaml"


# Summary
echo "CCDB_CONNECTION : $CCDB_CONNECTION"
echo "coatjava : $coatjava"
echo "yamlfile : $yamlfile"
echo "dec command  : $coatjava/coatjava/bin/decoder    -i $eviofile  -o $decfile"
echo "rec command  : $coatjava/coatjava/bin/recon-util -i $decfile   -o $recfile -y $yamlfile"

# Check CCDB
echo "> CCDB table check"
ccdb -r 22712 dump /calibration/alert/ahdc/time_offsets | grep "192 *0\.00"

# Uncomment to run command
$coatjava/coatjava/bin/decoder    -i "$eviofile"  -o "$decfile"
$coatjava/coatjava/bin/recon-util -i "$decfile"  -o "$recfile "-y "$yamlfile"
