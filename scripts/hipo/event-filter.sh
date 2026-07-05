#!/bin/bash

# filter elastics event for one file
# Make sure to modify the filtering criteria in io.github.ftouchte.filtering.ElasticEventFilter.java

# $1 is the input filename
# $2 is the output filename

home_dir=/w/hallb-scshelf2102/clas12/users/touchte

java -cp "$home_dir/coatjava/coatjava/lib/*:$home_dir/coatjava/coatjava/lib/clas/*:$home_dir/coatjava/coatjava/lib/services/*:$home_dir/coatjava/coatjava/lib/utils/*:$home_dir/amon/java-utils/target/classes" io.github.ftouchte.filtering.ElasticEventFilter $1 $2
