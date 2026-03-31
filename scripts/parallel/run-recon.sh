#!/bin/bash

if [ "$1" == "help" ]; then
    echo ""
    echo "---------"
    echo "  Usage"
    echo "---------"
    echo ""
    echo "  All fields are mandatory. One may have to set the path COATJAVA and the YAML file"
    echo ""
    echo -e "\033[1m   ./run-recon.sh \033[0m [folder]  [ncpu]"
    echo ""
    echo -e "  [folder]   folder containing the small files,  output of \033[1m ./run-reduce.sh \033[0m"
    echo "  [ncpu  ]   number of cpu, can match the number of files"
    echo ""
    exit 1
fi


folder=$1
ncpu=$2

coatjava_dir="/w/hallb-scshelf2102/clas12/users/touchte/coatjava-proj-only/"
#coatjava_dir=""
#coatjava_dir="/w/hallb-scshelf2102/clas12/users/touchte/coatjava/"
#coatjava_dir="/w/hallb-scshelf2102/clas12/users/touchte/clas12Tags/geometry_source/coatjava_src/"
#yaml_file="../ahdc_config.yaml"
yaml_file="../alert_clas12_config.yaml"

if [[ ! "$ncpu" =~ ^[0-9]+$ ]]; then
	echo "> $ncpu is not a number, it should correspond to the number of cpu"
	exit 1
fi

if [ ! -d "./$folder" ]; then
	echo "> This folder does not exists: $folder"
	exit 1
fi

#if [[ $ncpu -ge 32 ]]; then
#	echo "You attend to use more than 32 cpus: ncpu --> $ncpu"
#	echo "If you want to continue, you have to comment that verification in the script"
#	exit 1
#fi

#----------------------
# Start job
#----------------------
cd $folder
if [ ! -f "./$yaml_file" ]; then
	echo "> $yaml_file : file does not exist in this repository, it should be the yaml file"
	exit 1
fi
parallel --progress -j "$ncpu" "$coatjava_dir/coatjava/bin/recon-util -y $yaml_file -i {} -o rec-{/} > {/.}.log 2>&1" ::: file_*.hipo 

