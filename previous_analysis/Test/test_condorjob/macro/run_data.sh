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

this_script=$BASH_SOURCE
this_script=`readlink -f $this_script`
this_dir=`dirname $this_script`
echo rsyncing from $this_dir
echo running: $this_script $*

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
outDir="${!#}"

$list1 = $1
$list2 = $2
if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]
then
  cd $_CONDOR_SCRATCH_DIR
  rsync -av $this_dir/* .
else
  echo condor scratch NOT set
  exit -1
fi
  getinputfiles.pl --filelist $list1
  getinputfiles.pl --filelist $list2

# print the environment - needed for debugging
printenv

# this is how you run your Fun4All_G4_sPHENIX.C macro which was rsynced from your initial dir: 
root.exe -q -b Fun4All_FieldOnAllTrackersCalos.C\($nEvents,${inputFiles},\"${outDir}\"\)
#cp <output> /sphenix/tg/tg01/...

echo all done
