#!/usr/bin/bash

# creates file lists for the DST from known locations in lustre
# run number is the input argument

dir_streaming1=/sphenix/lustre01/sphnxpro/physics/slurp/streaming/physics/tpconlyrun_00049600_00049700
dir_streaming2=/sphenix/lustre01/sphnxpro/physics/slurp/streaming/physics/tpconlyrun_00049700_00049800
dir_calo1=/sphenix/lustre01/sphnxpro/physics/slurp/caloy2test/new_2024p006/run_00049600_00049700
dir_calo2=/sphenix/lustre01/sphnxpro/physics/slurp/caloy2test/new_2024p006/run_00049700_00049800
#dir_calo1=/sphenix/lustre01/sphnxpro/physics/slurp/caloy2calib/ana430_2024p007/run_00049600_00049700
#dir_calo2=/sphenix/lustre01/sphnxpro/physics/slurp/caloy2calib/ana430_2024p007/run_00049700_00049800

run=(49684 49685 49686 49687 49688 49689 49690 49691 49703 49704 49705 49706 49707 49708 49709 49710 49711 49712 49713 49714 49715 49716 49717 49718 49719 49720 49721 49722)

cd runList

for i in ${run[@]}
do
echo making run ${i} tpc list
if [ "$i" -lt "49700" ]; then
  echo $dir_streaming1
  ls -1 $dir_streaming1/*${i}* > run${i}_tpc.list
  if [ ! -s run${i}_tpc.list ]
  then
    echo run${i}_tpc.list empty, removing it
    rm run${i}_tpc.list
  fi
elif [ "$i" -ge "49700" ]; then
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
if [ "$i" -lt "49700" ]; then
  echo $dir_calo1
  ls -1 $dir_calo1/DST_CALO_*${i}* > run${i}_calo.list
  if [ ! -s run${i}_calo.list ]
  then
    echo run${i}_calo.list empty, removing it
    rm run${i}_calo.list
  fi
elif [ "$i" -ge "49700" ]; then
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
