#!/usr/bin/bash

# creates file lists for the DST from known locations in lustre
# run number is the input argument

dir_streaming=/sphenix/lustre01/sphnxpro/physics/slurp/streaming/physics/new_2024p002/run_00051300_00051400
dir_calo=/sphenix/lustre01/sphnxpro/physics/slurp/caloy2calib/ana430_2024p007/run_00051300_00051400

#full list
run=(51380 51381 51382 51383 51384 51385 51386 51387 51388 51389 51390 51391 51392 51394 51395 51396)

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

for i in ${run[@]}
do
if [ -f run${i}_calo.list ] && [ -f run${i}_tpc.list ]; then
    echo $i 2 >> ../run.list
    echo "Run ${i} have both tpc and calo dst"
else
    echo "Run ${i} don't have both tpc and calo dst"
fi
done

cd ..
