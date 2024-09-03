#!/usr/bin/bash

# creates file lists for the DST from known locations in lustre
# run number is the input argument

dir_streaming=/sphenix/lustre01/sphnxpro/physics/slurp/streaming/physics/tpconlyrun_00050900_00051000
dir_calo=/sphenix/lustre01/sphnxpro/physics/slurp/caloy2calib/ana430_2024p007/run_00050900_00051000

run=(50911 50913 50914 50915 50916 50918 50920 50921 50922 50924 50925 50927 50928 50929 50930 50932 50933 50934 50935 50936)

cd runList

for i in ${run[@]}
do
echo making run ${i} tpc list
ls -1 $dir_streaming/*${i}* > run${i}_tpc.list
if [ ! -s run${i}_tpc.list ]
then
  echo run${i}_tpc.list empty, removing it
  rm run${i}_tpc.list
fi
done

for i in ${run[@]}
do
echo making run ${i} calo list
ls -1 $dir_calo/DST_CALO_*${i}* > run${i}_calo.list
if [ ! -s run${i}_calo.list ]
then
  echo run${i}_calo.list empty, removing it
  rm run${i}_calo.list
fi
done

cd ..
