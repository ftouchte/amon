#!/bin/bash

outDir=$1

ccdb dump /calibration/alert/ahdc/time_to_distance_wire > $1/time2distance.txt

ccdb dump /geometry/alert/ahdc/layer_alignment > $1/layer_angles.txt

ccdb dump /geometry/alert/ahdc/wire_alignment > $1/wire_angles.txt
