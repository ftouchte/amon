#!/bin/bash

# --- working directory
workDir="/lustre24/expphy/volatile/clas12/touchte/alignment-cooking"

# --- Location of the inputs files
inputs="$workDir/decoded"
outputs="$workDir/reconstructed"

input_list="$workDir/input_list.txt"
output_list="$workDir/output_list.txt"

# Create the list of input files
find "$inputs" -regextype posix-extended -regex ".*/clas_[0-9]+\.evio\.[0-9]+\.hipo" | sort > "$input_list"

# Create the list of output files
while read -r file
do
    # Remove path for the file name
    filename=$(basename "$file")
    # extract the run number
    if [[ "$filename" =~ clas_([0-9]+)\. ]]; then
        runno="${BASH_REMATCH[1]}"
    fi
    # Print the output file in output_list.txt
    echo "${outputs}/$runno/rec_${filename}"
done < "$input_list" > "$output_list"

# Count number of files
N=$(wc -l < "$input_list")

echo "Number of files to be processed by SLURM:  $N"


# !!! The user should use the files created here (input_list.txt and output_list.txt and more... in the script run.sbatch)

# # Uncomment for testing
# for i in 0 1 2 3
# do
#   export SLURM_ARRAY_TASK_ID=$i
#   bash /w/hallb-scshelf2102/clas12/users/touchte/amon/scripts/slurm/run.sbatch
# done

# Uncomment to run command
#sbatch --array=0-$((N-1)) /w/hallb-scshelf2102/clas12/users/touchte/amon/scripts/slurm/run.sbatch