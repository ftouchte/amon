#!/bin/bash

runno=$(printf "%06d" "$1")
npart=$(printf "%05d" "$2")

# coatjava dir
coatjava="/w/hallb-scshelf2102/clas12/users/touchte"

#filename format
filename="clas_$runno.evio.$npart"

#eviofile="/cache/clas12/rg-l/data/clas_$runno/clas_$runno.evio.$npart"

# some paths
volatile="/lustre24/expphy/volatile/clas12/touchte/alignment-cooking"



decfile="$volatile/decoded/$runno/$filename.hipo"
recfile="$volatile/reconstructed/$runno/rec_$filename.hipo"

yamlfile="$volatile/config.yaml"

# echo "$coatjava/coatjava/bin/decoder    -i $eviofile     -o $decfile"
# $coatjava/coatjava/bin/decoder    -i $eviofile     -o $decfile


echo "$coatjava/coatjava/bin/recon-util -i $decfile  -o $recfile -y $yamlfile"
#$coatjava/coatjava/bin/recon-util -i "$decfile"  -o "$recfile" -y "$yamlfile"

