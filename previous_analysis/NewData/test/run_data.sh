#!/bin/bash

#source /sphenix/u/xyu3/.setup_sphenix.sh

source /opt/sphenix/core/bin/sphenix_setup.sh -n new  # setup sPHENIX environment in the singularity container shell. Note the shell is bash by default

# Additional commands for my local environment
export SPHENIX=/sphenix/u/xyu3
#export MYINSTALL=$SPHENIX/install
export MYINSTALL=/sphenix/u/xyu3/workarea/PhotonConv/module_test/install

# Setup MYINSTALL to local directory and run sPHENIX setup local script
# to adjust PATH, LD LIBRARY PATH, ROOT INCLUDE PATH, etc
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

echo "sPHENIX environment setup finished"

this_script=$BASH_SOURCE
this_script=`readlink -f $this_script`
this_dir=`dirname $this_script`
echo running: $this_script $*

root.exe -q -b Fun4All_FieldOnAllTrackersCalos.C
echo Script done
