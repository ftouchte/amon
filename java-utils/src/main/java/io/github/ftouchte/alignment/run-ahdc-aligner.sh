#!/bin/bash

# SQLITE_FILE="/work/clas12/users/touchte/ccdb_2026-05-24.sqlite"
SQLITE_FILE="/volatile/clas12/touchte/new-alignment/test_new_approach_layer/ccdb_2026-05-24.sqlite"

export CCDB_CONNECTION="sqlite:///$SQLITE_FILE"

if [[ ! -e $SQLITE_FILE ]]; then 
    echo "Sqlite file not found : $SQLITE_FILE"
    exit
else 
    echo "CCDB_CONNECTION = $CCDB_CONNECTION"
fi

home_dir=/w/hallb-scshelf2102/clas12/users/touchte

#java -cp "$home_dir/coatjava/coatjava/lib/*:$home_dir/amon/java-utils/target/classes" io.github.ftouchte.AhdcAlignmentAnalyser  "$@"
java -cp "$home_dir/coatjava/coatjava/lib/*:$home_dir/coatjava/coatjava/lib/clas/*:$home_dir/coatjava/coatjava/lib/services/*:$home_dir/coatjava/coatjava/lib/utils/*:$home_dir/amon/java-utils/target/classes" \
    io.github.ftouchte.alignment.AhdcAlignmentAnalyser "$@"
