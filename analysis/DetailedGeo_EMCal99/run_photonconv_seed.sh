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
InClusterDst=$2
InClusterPath=$3
InSeedDst=$4
InSeedPath=$5
InCaloDst=$6
InCaloPath=$7
OutPrefix=$8
OutPath=$9
Index=${10}
StepSize=${11}

root.exe -q -l Fun4All_TrackFitting_PhotonConv.C\($nEvents,\"${InClusterDst}\",\"${InClusterPath}\",\"${InSeedDst}\",\"${InSeedPath}\",\"${InCaloDst}\",\"${InCaloPath}\",\"${OutPrefix}\",\"${OutPath}\",$Index,$StepSize\)
echo Script done
