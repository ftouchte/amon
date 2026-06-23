#!/bin/bash

echo $CCDB_CONNECTION
echo "cmd : /calibration/alert/ahdc/time_to_distance_wire $1"

ccdb add -r - /calibration/alert/ahdc/time_to_distance_wire $1

echo "cmd : /calibration/alert/ahdc/time_to_distance_wire"

ccdb dump /calibration/alert/ahdc/time_to_distance_wire
