#!/bin/bash

runGroup="rgl"
coatjava="/work/clas12/users/touchte/coatjava/coatjava"
reconYaml="/volatile/clas12/touchte/alignment-cooking/config.yaml"
inputs="/mss/clas12/rg-l/data"
runs="22990"
outDir="/volatile/clas12/touchte/alignment-cooking"

# Create a minimaliste config file
cat <<EOF > cfg.json
{
  "runGroup" : "$runGroup",
  "coatjava" : "$coatjava",
  "mergeSize" : 5,
  "phaseSize" : 5000,
  "reconYaml" : "$reconYaml",
  "inputs" : "$inputs",
  "runs" : "$runs",
  "clara": "/scigroup/cvmfs/hallb/clas12/sw/noarch/clara/5.0.2_13.2.0"
}

EOF


# execute command
 clas12-workflow --config cfg.json --tag ahdc_wire_position --model decrec \
 --runs 22990-22994 \
 --outDir $outDir  \
 --reconYaml $reconYaml \
 --submit




