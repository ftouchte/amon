#!/bin/bash

ncpu=4
coatjava_dir="/w/hallb-scshelf2102/clas12/users/touchte/coatjava"

#while IFS= read -r line; do
#    echo "$line"
#    # Read input file
#    file="/cache/clas12/rg-l/data/clas_022712/$line"
#    printf "\033[1;32m INPUT\033[0m $file\n"
#
#    # Decode input file
#    decfile="dec_$line.hipo"
#    echo "$coatjava_dir/coatjava/bin/decoder -i $file -o $decfile"
#
#    # Run recon on this file
#    recfile="rec_$line.hipo"
#    echo "$coatjava_dir/coatjava/bin/recon-util -i $decfile -o $recfile"
#
#done < eviofiles.txt


#parallel --dry-run --progress -j 1 '
parallel --progress -j $ncpu '
    file="/cache/clas12/rg-l/data/clas_022712/{}"

    printf "\033[1;32mINPUT\033[0m %s\n" "$file"

    decfile="dec/dec_{/}.hipo"
    recfile="rec/rec_{/}.hipo"

    declog="log/dec_{/}.log"
    reclog="log/rec_{/}.log"
    
    '"$coatjava_dir"'/coatjava/bin/decoder -i "$file" -o "$decfile" > "$declog" 2>&1

    '"$coatjava_dir"'/coatjava/bin/recon-util -y /w/hallb-scshelf2102/clas12/users/touchte/data/simu/alert_clas12_config.yaml -i "$decfile" -o "$recfile" > "reclog" 2>&1
' :::: ../eviofiles.txt
