#!/bin/bash

# filter elastics events for all files in a given directory

path=$1

outdir="$path/elastic-filtered"

if [[ ! -d $outdir ]]; then 
    mkdir "$outdir"
fi

for file in "$path"/rec*
#for file in v120/rec_clas_022712.evio.0000[0-3].hipo
do
    echo -e "\033[1;32m > Filter elastics for : \033[0m $file"
    filename=$(basename "$file") # remove the path
    #name="${filename%.*}" # remove the extension
    /w/hallb-scshelf2102/clas12/users/touchte/amon/scripts/hipo/event-filter.sh "$file" "$outdir/$filename"
done
