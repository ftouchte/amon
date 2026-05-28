#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --job-name=run_ahdc_wire_alignment
#SBATCH --output=/farm_out/%u/%x-%A-%a.out
#SBATCH --error=/farm_out/%u/%x-%A-%a.err
#SBATCH --time=30:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --mem=40G

/w/hallb-scshelf2102/clas12/users/touchte/amon/scripts/hipo/run-ahdc-aligner.sh -i /volatile/clas12/touchte/new-translation-table/reconstructed/022712/elastic-filtered/merged/rec_clas_022712.evio.0-834.hipo -o /lustre24/expphy/volatile/clas12/touchte/new-alignment/wire_alignment_following -ncpu 60