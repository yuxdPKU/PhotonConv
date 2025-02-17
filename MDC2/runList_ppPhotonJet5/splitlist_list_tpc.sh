#!/bin/bash

input_file=$1
DST_per_job=$2

if [ ! -f "$input_file" ]; then
  echo "File not found: $input_file"
  exit 1
fi


while IFS= read -r runtype; do
  tracklist=${runtype}.list
  echo reading run ${tracklist}

  if [[ ! -f "$tracklist" ]]; then
    echo "'$tracklist' does not exist. Skipping..."
    continue
  fi

  line_number=0
  file_count=0

  while IFS= read -r line; do
    if (( line_number % $DST_per_job == 0 )); then
      if (( line_number != 0 )); then
        file_count=$((file_count + 1))
      fi
      output_file_dst=$(printf "trackrunlist/%s_%05d.list" "${runtype}" "$file_count")
      if [ -f "$output_file_dst" ]; then
        rm "$output_file_dst"
      fi
    fi

    echo "$line" >> "$output_file_dst"

    line_number=$((line_number + 1))
  done < "$tracklist"
  file_count=$((file_count + 1))
 
done < "$input_file"

