#!/bin/bash

input_file=$1
output_file=$2

if [ ! -f "$input_file" ]; then
  echo "File not found: $input_file"
  exit 1
fi

while IFS= read -r runnumber; do
  echo reading run $runnumber
  trkrtrackslist=dstlist/dst_trkr_tracks_run2pp-000${runnumber}.list
  trkrclusterlist=dstlist/dst_trkr_cluster_run2pp-000${runnumber}.list
  calolist=dstlist/dst_calo_run2pp-000${runnumber}.list
  #streaminglist=dst_streaming_event_run2pp-000${runnumber}.list

  if [[ ! -f "$trkrtrackslist" || ! -f "$trkrclusterlist" || ! -f "$calolist" ]]; then
    echo "Either '$trkrtrackslist' or '$trkrclusterlist' or '$calolist' does not exist. Skipping..."
    continue
  fi

  echo ${runnumber} >> ${output_file}
 
done < "$input_file"

cd ../
