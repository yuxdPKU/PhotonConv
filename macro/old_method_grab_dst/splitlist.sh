#!/bin/bash

#input_file=$1
#run_number=$2

run_number=$1
input_file=run${run_number}_tpc.list

if [ ! -f "$input_file" ]; then
  echo "File not found: $input_file"
  exit 1
fi

line_number=0

while IFS= read -r line; do
  output_file=$(printf "trackrunlist/run%d_%04d.txt" "$run_number" "$line_number")
  
  echo "$line" > "$output_file"
  
  line_number=$((line_number + 1))
done < "$input_file"

