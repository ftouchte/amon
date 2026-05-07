#!/bin/bash

# filter elastics event for one file
# Make sure to modify the filtering criteria in io.github.ftouchte.filtering.ElasticEventFilter.java

java -cp /w/hallb-scshelf2102/clas12/users/touchte/coatjava/coatjava/lib/clas/coat-libs-13.7.1-SNAPSHOT.jar:/w/hallb-scshelf2102/clas12/users/touchte/amon/java-utils/target/classes io.github.ftouchte.filtering.ElasticEventFilter $1 $2
