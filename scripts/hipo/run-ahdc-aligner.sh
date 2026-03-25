#!/bin/bash

home_dir=/w/hallb-scshelf2102/clas12/users/touchte

#java -cp "$home_dir/coatjava/coatjava/lib/*:$home_dir/amon/java-utils/target/classes" io.github.ftouchte.AhdcAlignmentAnalyser  "$@"
java -cp "$home_dir/coatjava/coatjava/lib/*:$home_dir/coatjava/coatjava/lib/clas/*:$home_dir/coatjava/coatjava/lib/services/*:$home_dir/coatjava/coatjava/lib/utils/*:$home_dir/amon/java-utils/target/classes" \
    io.github.ftouchte.AhdcAlignmentAnalyser "$@"
