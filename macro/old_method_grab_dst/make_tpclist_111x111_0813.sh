#!/usr/bin/bash

# creates file lists for the DST from known locations in lustre
# run number is the input argument

dir_streaming=/sphenix/lustre01/sphnxpro/physics/slurp/streaming/physics/new_2024p002/run_00051200_00051300
dir_calo=/sphenix/lustre01/sphnxpro/physics/slurp/caloy2calib/ana430_2024p007/run_00051200_00051300

run=(51240 51241 51242 51243 51244 51245 51246 51247 51248 51249 51250 51252 51253 51254 51255 51276 51277 51279 51281 51282)

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
