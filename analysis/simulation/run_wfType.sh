#!/bin/bash

./wfType.exe "/home/touchte-codjo/Desktop/hipofiles/wfType/redecoded-raw-D2-run-23003.hipo"
cp ./output/wfType_study.root ./output/wfType_study-raw.root
./wfType.exe "/home/touchte-codjo/Desktop/hipofiles/wfType/redecoded-elastics-proton.hipo"
cp ./output/wfType_study.root ./output/wfType_study-proton.root
./wfType.exe "/home/touchte-codjo/Desktop/hipofiles/wfType/redecoded-elastics-deuteron.hipo"
cp ./output/wfType_study.root ./output/wfType_study-deuteron.root

ls -l ./output/wfType_study*
