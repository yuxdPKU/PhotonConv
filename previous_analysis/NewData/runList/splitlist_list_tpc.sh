#!/bin/bash

input_file=$1
output_file=$2
DST_per_job=$3

if [ ! -f "$input_file" ]; then
  echo "File not found: $input_file"
  exit 1
fi


while IFS= read -r runnumber; do
  echo reading run $runnumber
  streaminglist=dst_streaming_event_run2pp-000${runnumber}.list
  #streaminglist=dst_tpc_event_run2pp-000${runnumber}.list
  calolist=dst_calo_run2pp-000${runnumber}.list

  if [[ ! -f "$streaminglist" || ! -f "$calolist" ]]; then
    echo "Either '$streaminglist' or '$calolist' does not exist. Skipping..."
    continue
  fi

  line_number=0
  file_count=0

  while IFS= read -r line; do
    if (( line_number % $DST_per_job == 0 )); then
      if (( line_number != 0 )); then
        file_count=$((file_count + 1))
      fi
      output_file_dst=$(printf "trackrunlist/run%d_%04d.list" "$runnumber" "$file_count")
      if [ -f "$output_file_dst" ]; then
        rm "$output_file_dst"
      fi
    fi

    echo "$line" >> "$output_file_dst"

    line_number=$((line_number + 1))
  done < "$streaminglist"
  file_count=$((file_count + 1))
  echo ${runnumber} ${file_count} >> ${output_file}
 
done < "$input_file"

