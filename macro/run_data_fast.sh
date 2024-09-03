#!/bin/bash

source /sphenix/u/xyu3/.setup_sphenix.sh

nEvents=-1
#nEvents=100

count=1
total=$#

inputFiles="{"
for fileList in $@
do
  if [ $count -eq $total ]; then
    break
  fi
  inputFiles+="\"${fileList}\","
  count=$((count + 1))
done
inputFiles=${inputFiles::-1}
inputFiles+="}"
DrawEvtDisplay="${!#}"
echo running: run_data_fast.sh $*
#valgrind --num-callers=30 --leak-check=full --suppressions=$ROOTSYS/etc/valgrind-root.supp root.exe -q -b Fun4All_Template.C\(${inputFiles},$nEvents\)
#valgrind --num-callers=30 --leak-check=full --show-leak-kinds=all --track-origins=yes --suppressions=$ROOTSYS/root.supp root.exe -q -b Fun4All_Template.C\(${inputFiles},$nEvents\)
root.exe -q -b Fun4All_TrackAnalysis.C\($nEvents,${inputFiles},true,${DrawEvtDisplay}\)
echo Script done
