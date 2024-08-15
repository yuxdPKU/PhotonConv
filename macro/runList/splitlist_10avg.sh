#!/bin/bash

run_number=$1
input_file=run${run_number}_tpc.list

if [ ! -f "$input_file" ]; then
  echo "File not found: $input_file"
  exit 1
fi

line_number=0
file_count=0

while IFS= read -r line; do
  if (( line_number % 10 == 0 )); then
    if (( line_number != 0 )); then
      file_count=$((file_count + 1))
    fi
    output_file=$(printf "trackrunlist/run%d_%04d.txt" "$run_number" "$file_count")
  fi

  echo "$line" >> "$output_file"

  line_number=$((line_number + 1))
done < "$input_file"

