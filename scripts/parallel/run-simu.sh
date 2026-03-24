#!/bin/bash

if [ "$1" == "help" ]; then
    echo ""
    echo "---------"
    echo "  Usage"
    echo "---------"
    echo ""
    echo "  All fields are mandatory. The output file can be merged later using hipo-utils"
    echo ""
    echo -e "\033[1m   ./run-simu.sh \033[0m [folder]  [ncpu]  [runno]  [outdir]"
    echo ""
    echo "  [folder]   folder containing the lund files"
    echo "  [ncpu  ]   number of cpu, can match the number of lund files"
    echo "  [runno ]   run number"
    echo "  [outdir]   output directory"
    echo ""
    exit 1
fi

# this folder contains lund files
folder=$1
ncpu=$2
runno=$3
outdir=$4

gemc_dir="/w/hallb-scshelf2102/clas12/users/touchte/clas12Tags/source"
gcard_file="../alert-lund.gcard"

if [ ! -d "$folder" ]; then
	echo "> This folder does not exists: $folder"
	exit 1
fi

if [[ ! "$ncpu" =~ ^[0-9]+$ ]]; then
	echo "> $ncpu is not a number, it should correspond to the number of cpu"
	exit 1
fi

if [[ ! "$runno" =~ ^[0-9]+$ ]]; then
	echo "> $runno is not a number, it should correspond to the run number"
	exit 1
fi

if [ ! -d "./$outdir" ]; then
	echo "> This folder does not exists: $folder"
    mkdir $outdir
    echo "New folder created: $outdir"
fi
outdir=$(realpath $outdir)


#if [[ $ncpu -ge 32 ]]; then
#	echo "You attend to use more than 32 cpus: ncpu --> $ncpu"
#	echo "If you want to continue, you have to comment that verification in the script"
#	exit 1
#fi

#----------------------
# Start job
#----------------------
cd $folder
if [ ! -f "./$gcard_file" ]; then
	echo "> $gcard_file : file does not exist in this repository, it should be the gcard file"
	exit 1
fi
#parallel --dry-run --progress -j "$ncpu" "$gemc_dir/gemc $gcard_file -USE_GUI=0 -RUNNO=$runno -N=$nevents -INPUT_GEN_FILE='LUND, {/}' -OUTPUT='hipo, simu-{/.}.hipo' > {/.}.log 2>&1" ::: lund_file_*.dat 
parallel --progress -j "$ncpu" "$gemc_dir/gemc $gcard_file -USE_GUI=0 -RUNNO=$runno -INPUT_GEN_FILE='LUND, {/}' -OUTPUT='hipo, $outdir/file-{/.}.hipo' > $outdir/{/.}.log 2>&1" ::: lund_file_*.dat

