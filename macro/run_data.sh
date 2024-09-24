#!/bin/bash

#source /sphenix/u/xyu3/.setup_sphenix.sh

source /opt/sphenix/core/bin/sphenix_setup.sh -n new  # setup sPHENIX environment in the singularity container shell. Note the shell is bash by default

# Additional commands for my local environment
export SPHENIX=/sphenix/u/xyu3
export MYINSTALL=$SPHENIX/install

# Setup MYINSTALL to local directory and run sPHENIX setup local script
# to adjust PATH, LD LIBRARY PATH, ROOT INCLUDE PATH, etc
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

echo "sPHENIX environment setup finished"

nEvents=$1

#skip the first parameter
shift

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
echo running: run_data.sh $*
#valgrind --num-callers=30 --leak-check=full --suppressions=$ROOTSYS/etc/valgrind-root.supp root.exe -q -b Fun4All_Template.C\(${inputFiles},$nEvents\)
#valgrind --num-callers=30 --leak-check=full --show-leak-kinds=all --track-origins=yes --suppressions=$ROOTSYS/root.supp root.exe -q -b Fun4All_Template.C\(${inputFiles},$nEvents\)
root.exe -q -b Fun4All_FieldOnAllTrackersCalos.C\($nEvents,${inputFiles},true,${DrawEvtDisplay}\)
echo Script done
