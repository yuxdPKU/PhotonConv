#!/usr/bin/bash

# creates file lists for the DST from known locations in lustre
# run number is the input argument

dir_streaming1=/sphenix/lustre01/sphnxpro/physics/slurp/streaming/physics/tpconlyrun_00051000_00051100
dir_streaming2=/sphenix/lustre01/sphnxpro/physics/slurp/streaming/physics/tpconlyrun_00051100_00051200
dir_calo1=/sphenix/lustre01/sphnxpro/physics/slurp/caloy2calib/ana430_2024p007/run_00051000_00051100
dir_calo2=/sphenix/lustre01/sphnxpro/physics/slurp/caloy2calib/ana430_2024p007/run_00051100_00051200

run=(51099 51100 51101 51102 51103 51104 51105 51119 51106 51107 51109)

cd runList

for i in ${run[@]}
do
echo making run ${i} tpc list
if [ "$i" -lt "51100" ]; then
  echo $dir_streaming1
  ls -1 $dir_streaming1/*${i}* > run${i}_tpc.list
  if [ ! -s run${i}_tpc.list ]
  then
    echo run${i}_tpc.list empty, removing it
    rm run${i}_tpc.list
  fi
elif [ "$i" -ge "51100" ]; then
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
if [ "$i" -lt "51100" ]; then
  echo $dir_calo1
  ls -1 $dir_calo1/DST_CALO_*${i}* > run${i}_calo.list
  if [ ! -s run${i}_calo.list ]
  then
    echo run${i}_calo.list empty, removing it
    rm run${i}_calo.list
  fi
elif [ "$i" -ge "51100" ]; then
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
