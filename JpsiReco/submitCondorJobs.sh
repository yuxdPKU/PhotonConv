#!/bin/bash

submissionFile=condor-data.job

#runs="53877"
#runs="53046 53079 53081 53195 53494 53513 53517 53530 53531 53532 53571 53578 53580 53581 53590 53630 53631 53632 53652 53687 53716 53742 53744 53756" #green runs, both good Tracking and calo runs
runs="53018 53194 53196 53197 53586 53587 53741 53743 53783 53871 53876 53877" #yellow runs
#runs="53080 53738 53739 53879" #red runs

for run in $runs
do

  lower_bound=$(( (${run} / 100) * 100 ))
  upper_bound=$(( lower_bound + 100 ))

  formatted_lower=$(printf "%08d" "$lower_bound")
  formatted_upper=$(printf "%08d" "$upper_bound")

  sed -i -e "s/^RunNumber[[:space:]]*=[[:space:]]*[0-9]*/RunNumber      = $run/" \
       -e "s/^RunRange[[:space:]]*=[[:space:]]*run_[0-9]*_[0-9]*/RunRange       = run_${formatted_lower}_${formatted_upper}/" ${submissionFile}

  ListFile=filelist/dst_sync_trkr_calo_run2pp-$(printf "%08d" "$run").list
  if [ -f ${ListFile} ]; then
    echo "New list file: ${ListFile}"
    condor_submit ${submissionFile}
  else
    echo "List lile ${ListFile} is missing"
  fi

done
