#!/bin/bash

dossier=$1
version=$2
nature=$3


if [ ! -d "./$dossier" ]; then
	echo "> This folder does not exist: $dossier"
	exit 1
fi

if [[ ! "$version" =~ ^[0-9]+$ ]]; then
        echo "> $version is not a number, it should correspond to the version number"
        exit 1
fi
#----------------------
# Start job
#----------------------

cd $dossier

if [ $nature = "data" ]; then  
    hipo-utils -merge -o "rec-data-r22712-v$version.hipo" rec-file*
    echo ">> Path of the new file:"
    ls $PWD/rec-data*
else
    hipo-utils -merge -o "rec-simu-deuteron-v$version.hipo" rec-file*
    echo ">> Path of the new file:"
    ls $PWD/rec-simu*
fi


