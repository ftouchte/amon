#!/bin/bash

if [ "$1" == "help" ]; then
    echo ""
    echo "---------"
    echo "  Usage"
    echo "---------"
    echo ""
    echo "  All fields are mandatory. Use \`hipo-utils -split ...\`. Example of output name: file_1.hipo. To be used with run recon."
    echo ""
    echo -e "\033[1m   ./run-reduce.sh \033[0m [filename]  [nevents]  [outdir]"
    echo ""
    echo "  [filename]   file to be reduced in small parts"
    echo "  [nevents ]   number of events in each files"
    echo "  [outdir  ]   output directory"
    echo ""
    exit 1
fi

filename=$1
nevents=$2
dossier=$3
#ncpu=$3

#hipo-utils -info $filename | awk '/number of  events/ {print $NF}'

#set -e


if [ ! -f "./$filename" ]; then
	echo "> $filename does not exist in this repository"
	exit 1
fi

if [[ ! "$nevents" =~ ^[0-9]+$ ]]; then
	echo "> $nevents is not a number"
	exit 1
fi

if [ -d "./$dossier" ]; then
	echo "> This folder already exists: $dossier"
	exit 1
fi

#----------------------
# Start job
#----------------------

echo "> Processing $filename..."
echo "> Run: hipo-utils -split -n $nevents $filename"
	hipo-utils -split -n $nevents $filename
echo "> Folder created: $dossier"
	num="0";
	files=(${filename}_*)
	mkdir $dossier
	cd $dossier
	#ls ${filename}_*
	for file in "${files[@]}"; do
		((num++))
		if [ -e "../$file" ]; then
			mv "../$file" "file_$num.hipo"
			#echo "> new file file_$num.hipo"
		fi
	done
	echo "Source file: ../$filename" > readme.txt
	echo "Number of portions: $num" >> readme.txt
echo "> All files moved in: $dossier"
echo "Number of portions: $num"

