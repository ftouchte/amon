#!/bin/bash

coatjava="/w/hallb-scshelf2102/clas12/users/touchte/coatjava-kf-fit"

runno=$1
npart=$(printf "%05d" $2)

filename="clas_0$runno.evio.$npart"

eviofile="/cache/clas12/rg-l/data/clas_0$runno/clas_0$runno.evio.$npart"

volatiledir="/lustre24/expphy/volatile/clas12/touchte"

decfile="$volatiledir/22267/dec/dec_$filename.hipo"

recfile="$volatiledir/22267/rec/rec_$filename.hipo"

yamlfile="$volatiledir/22267/config.yaml"

    echo -e "\033[1;32m> Decoding \033[0m"
    echo "$coatjava/coatjava/bin/decoder    -i $eviofile     -o $decfile"
$coatjava/coatjava/bin/decoder    -i $eviofile     -o $decfile

    echo -e "\033[1;32m> Reconstruction with KF fit \033[0m"
    echo "$coatjava/coatjava/bin/recon-util -i $decfile  -o $recfile -y $yamlfile"
$coatjava/coatjava/bin/recon-util -i $decfile  -o $recfile -y $yamlfile

