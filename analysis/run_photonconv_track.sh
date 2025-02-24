#!/bin/bash

source /opt/sphenix/core/bin/sphenix_setup.sh -n new  # setup sPHENIX environment in the singularity container shell. Note the shell is bash by default

# Additional commands for my local environment
export SPHENIX=/sphenix/u/xyu3
export MYINSTALL=$SPHENIX/install

# Setup MYINSTALL to local directory and run sPHENIX setup local script
# to adjust PATH, LD LIBRARY PATH, ROOT INCLUDE PATH, etc
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

echo "sPHENIX environment setup finished"

this_script=$BASH_SOURCE
this_script=`readlink -f $this_script`
this_dir=`dirname $this_script`
echo running: $this_script $*

nEvents=$1
InTrackDst=$2
InTrackPath=$3
InCaloDst=$4
InCaloPath=$5
OutPrefix=$6
OutPath=$7
RunNumber=$8
Segment=$9
Index=${10}
StepSize=${11}

root.exe -q -b Fun4All_TrackAnalysis_PhotonConv.C\($nEvents,\"${InTrackDst}\",\"${InTrackPath}\",\"${InCaloDst}\",\"${InCaloPath}\",\"${OutPrefix}\",\"${OutPath}\",$RunNumber,$Segment,$Index,$StepSize\)
echo Script done
