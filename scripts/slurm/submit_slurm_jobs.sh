#!/bin/bash

##########################################################################
# Usage: better to read the code
#      !!!! Need to modify the complementary code run.sbatch acoordingly
##########################################################################

# --- working directory
workDir="/lustre24/expphy/volatile/clas12/touchte/new-translation-table" # should macth the one in run.sbatch
cache="/cache/clas12/rg-l/data"

# --- Location of the reconstructed files
decDir="$workDir/decoded"
recDir="$workDir/reconstructed"

if [[ ! -d "$decDir" ]]; then
    mkdir "$decDir"
fi

if [[ ! -d "$recDir" ]]; then
    mkdir "$recDir"
fi

evio_list="$workDir/evio_list.txt"
dec_list="$workDir/dec_list.txt"
rec_list="$workDir/rec_list.txt"

# --- Create the list of input (here evio) files
# find $cache/clas_02299[0124] -regextype posix-extended -regex ".*/clas_[0-9]+\.evio\.[0-9]+\.hipo" | sort > "$dec_list" # case of hipo file
find $cache/clas_02299[0124] -regextype posix-extended -regex ".*/clas_[0-9]+\.evio\.[0-9]+" | sort > "$evio_list"

# --- Create the list of dec files
while read -r file
do
    # Remove path for the file name
    filename=$(basename "$file")
    # extract the run number
    if [[ "$filename" =~ clas_([0-9]+)\. ]]; then
        runno="${BASH_REMATCH[1]}"
    fi
    # create the repository for this run number if it does not exist
    if [[ ! -d "$decDir/$runno" ]]; then 
        mkdir "$decDir/$runno"
    fi
    # Print the name of the decoded file
    echo "${decDir}/$runno/${filename}.hipo"
done < "$evio_list" > "$dec_list"


# --- Create the list of rec files
while read -r file
do
    # Remove path for the file name
    filename=$(basename "$file")
    # extract the run number
    if [[ "$filename" =~ clas_([0-9]+)\. ]]; then
        runno="${BASH_REMATCH[1]}"
    fi
    # create the repository for this run number if it does not exist
    if [[ ! -d "$recDir/$runno" ]]; then 
        mkdir "$recDir/$runno"
    fi
    # Print the name of the reconstructed file
    echo "${recDir}/$runno/rec_${filename}"
done < "$dec_list" > "$rec_list"

# Count number of files
N=$(wc -l < "$evio_list")

echo "Number of files to be processed by SLURM:  $N"


# !!! The user should use the files created here (input_list.txt and output_list.txt and more... in the script run.sbatch)

# Uncomment for testing : deterministics
# for i in 0 1 2 3
# do
#   export SLURM_ARRAY_TASK_ID=$i
#   bash /w/hallb-scshelf2102/clas12/users/touchte/amon/scripts/slurm/run.sbatch
# done

# Uncomment for testing : random
# for _ in 0 1 2 3
# do
#     export SLURM_ARRAY_TASK_ID=$((RANDOM % N))
#     bash /w/hallb-scshelf2102/clas12/users/touchte/amon/scripts/slurm/run.sbatch
# done



# Uncomment to run command
sbatch --array=0-$((N-1)) /w/hallb-scshelf2102/clas12/users/touchte/amon/scripts/slurm/run.sbatch