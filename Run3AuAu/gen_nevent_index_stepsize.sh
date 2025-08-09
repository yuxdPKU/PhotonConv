#!/bin/bash

output_file="nevent_index_stepsize.txt"

> "$output_file"

for ((i=0; i<100; i++))
do
    echo "100 $i 100" >> "$output_file"
done

echo done
