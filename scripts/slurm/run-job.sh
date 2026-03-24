#!/bin/bash

coatjavav119="/w/hallb-scshelf2102/clas12/users/touchte/coatjava-proj-only"
coatjavav120="/w/hallb-scshelf2102/clas12/users/touchte/coatjava-kf-fit"

runno=$1
npart=$(printf "%05d" $2)

filename="clas_0$runno.evio.$npart"
eviofile="/cache/clas12/rg-l/data/clas_0$runno/clas_0$runno.evio.$npart"

volatiledir="/lustre24/expphy/volatile/clas12/touchte"

decfile="$volatiledir/kalman-filter/v119/dec_$filename.hipo"

recfilev119="$volatiledir/kalman-filter/v119/rec_$filename.hipo"
recfilev120="$volatiledir/kalman-filter/v120/rec_$filename.hipo"

yamlfilev119="$volatiledir/kalman-filter/v119/config.yaml"
yamlfilev120="$volatiledir/kalman-filter/v120/config.yaml"

echo -e "\033[1;32m> Decoding \033[0m"
#echo "$coatjavav119/coatjava/bin/decoder    -i $eviofile     -o $decfile"
#$coatjavav119/coatjava/bin/decoder    -i $eviofile     -o $decfile -n 10000
$coatjavav119/coatjava/bin/decoder    -i $eviofile     -o $decfile
echo -e "\033[1;32m> Reconstruction with KF propagation only \033[0m"
#echo "$coatjavav119/coatjava/bin/recon-util -i $decfile      -o $recfilev119 -y $yamlfilev119"
#$coatjavav119/coatjava/bin/recon-util -i $decfile      -o $recfilev119 -y $yamlfilev119 -n 10000
$coatjavav119/coatjava/bin/recon-util -i $decfile      -o $recfilev119 -y $yamlfilev119
echo -e "\033[1;32m> Reconstruction with KF fit \033[0m"
#echo "$coatjavav120/coatjava/bin/recon-util -i $recfilev119  -o $recfilev120 -y $yamlfilev120"
#$coatjavav120/coatjava/bin/recon-util -i $recfilev119  -o $recfilev120 -y $yamlfilev120 -n 10000
$coatjavav120/coatjava/bin/recon-util -i $recfilev119  -o $recfilev120 -y $yamlfilev120

