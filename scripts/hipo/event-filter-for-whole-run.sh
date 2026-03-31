#!/bin/bash

path=/lustre24/expphy/volatile/clas12/touchte/kalman-filter/v120
outdir=/lustre24/expphy/volatile/clas12/touchte/kalman-filter/elastic-filtered/v120

for file in $path/rec*
#for file in v120/rec_clas_022712.evio.0000[0-3].hipo
do
    echo -e "\033[1;32m > Filter elastics for : \033[0m $file"
    filename=$(basename "$file") # remove the path
    #name="${filename%.*}" # remove the extension
    /w/hallb-scshelf2102/clas12/users/touchte/amon/scripts/hipo/event-filter.sh $file "$outdir/$filename"
done
