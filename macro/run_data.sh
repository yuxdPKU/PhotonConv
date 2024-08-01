#!/bin/bash

source /sphenix/u/xyu3/.setup_sphenix.sh

nEvents=-1

inputFiles="{"
for fileList in $@
do
  inputFiles+="\"${fileList}\","
done
inputFiles=${inputFiles::-1}
inputFiles+="}"
echo running: run_data.sh $*
#valgrind --num-callers=30 --leak-check=full --suppressions=$ROOTSYS/etc/valgrind-root.supp root.exe -q -b Fun4All_Template.C\(${inputFiles},$nEvents\)
root.exe -q -b Fun4All_Template.C\(${inputFiles},$nEvents\)
echo Script done
