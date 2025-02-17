#!/bin/bash

directory="/sphenix/u/xyu3/hftg01/TPC_only_track_DST/Reconstructed"

if [ ! -d "$directory" ]; then
    echo "don't exist"
    exit 1
fi

for folder in "$directory"/*; do
    if [ -d "$folder" ]; then
        folder_name=$(basename "$folder")
        output_file="${folder_name}.list"
        
        > $output_file
        
#        for file in "$folder"/*; do
#            if [ -f "$file" ]; then
#                echo "$folder"/"$(basename $file)" >> $output_file
#            fi
#        done
        for file in $(ls "$folder" | sort -t'-' -k2 -n); do
            if [ -f "$folder/$file" ]; then
                echo "$folder"/"$file" >> $output_file
            fi
        done
    fi
done

echo "done"
