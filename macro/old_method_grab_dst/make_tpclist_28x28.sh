#!/usr/bin/bash

# creates file lists for the DST from known locations in lustre
# run number is the input argument

#dir_streaming1=/sphenix/lustre01/sphnxpro/physics/slurp/streaming/physics/tpconlyrun_00049900_00050000
#dir_streaming2=/sphenix/lustre01/sphnxpro/physics/slurp/streaming/physics/tpconlyrun_00050000_00051000
dir_streaming1=/sphenix/lustre01/sphnxpro/physics/slurp/streaming/fast/run_00049900_00050000
dir_streaming2=/sphenix/lustre01/sphnxpro/physics/slurp/streaming/fast/run_00050000_00050100
dir_calo1=/sphenix/lustre01/sphnxpro/physics/slurp/caloy2calib/ana430_2024p007/run_00049900_00050000
dir_calo2=/sphenix/lustre01/sphnxpro/physics/slurp/caloy2calib/ana430_2024p007/run_00050000_00050100

run=(49997 49998 49999 50000 50001 50002 50003 50004 50005 50006 50007 50008 50009 50010 50011 50012 50013 50014 50015 50016 50017 50018 50019 50020 50021 50022 50023 50024 50025)

cd runList

for i in ${run[@]}
do
echo making run ${i} tpc list
if [ "$i" -lt "50000" ]; then
  echo $dir_streaming1
  ls -1 $dir_streaming1/*${i}* > run${i}_tpc.list
  if [ ! -s run${i}_tpc.list ]
  then
    echo run${i}_tpc.list empty, removing it
    rm run${i}_tpc.list
  fi
elif [ "$i" -ge "50000" ]; then
  echo $dir_streaming2
  ls -1 $dir_streaming2/*${i}* > run${i}_tpc.list
  if [ ! -s run${i}_tpc.list ]
  then
    echo run${i}_tpc.list empty, removing it
    rm run${i}_tpc.list
  fi
fi
done

for i in ${run[@]}
do
echo making run ${i} calo list
if [ "$i" -lt "50000" ]; then
  echo $dir_calo1
  ls -1 $dir_calo1/DST_CALO_*${i}* > run${i}_calo.list
  if [ ! -s run${i}_calo.list ]
  then
    echo run${i}_calo.list empty, removing it
    rm run${i}_calo.list
  fi
elif [ "$i" -ge "50000" ]; then
  echo $dir_calo2
  ls -1 $dir_calo2/DST_CALO_*${i}* > run${i}_calo.list
  if [ ! -s run${i}_calo.list ]
  then
    echo run${i}_calo.list empty, removing it
    rm run${i}_calo.list
  fi
fi
done

cd ..
