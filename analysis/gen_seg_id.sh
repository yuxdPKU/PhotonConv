#!/bin/bash

#53877 - 400khz
#53876 - 430khz
#53756 - 380khz
#53744 - 300khz
#53630 - 550khz
#53534 - 250khz

get_closest_numbers() {
    number=$1
    lower=$(( (number / 100) * 100 ))
    upper=$(( lower + 100 ))
    echo "$lower $upper"
}

#runs=(53744)
runs=(53756)

path=/sphenix/u/xyu3/workarea/KsReco_sPHENIX/production/filelist

NumEvtPerDst=10000
NumEvtPerJob=2000
NumJobPerDst=$((${NumEvtPerDst} / ${NumEvtPerJob}))

for ((k=0; k<${#runs[@]}; k++))
do
  out=${runs[$k]}_seg_id.txt
  > ${out}

  nseg=`ls ${path}/${runs[$k]}/rawhit_*.list | wc -l`

  echo run ${runs[$k]} : ${nseg} dsts, $((${nseg} * ${NumJobPerDst})) jobs

  for ((i=0; i<${nseg}; i++))
  do
    calosegment=$((i / 10))
    for ((j=0; j<${NumJobPerDst}; j++))
    do
      echo "$i $j $calosegment" >> ${out}
    done
  done

done
