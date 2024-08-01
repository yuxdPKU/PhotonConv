#!/bin/bash
FILES=""
MERGE="final_46730_ana.root"
for f in TrackCalo_*_ana.root
do
	echo "Processing $f"
  FILES+="$f "
done
hadd -f $MERGE $FILES

FILES=""
MERGE="final_46730_kfp.root"
for f in TrackCalo_*_kfp.root
do
	echo "Processing $f"
  FILES+="$f "
done
hadd -f $MERGE $FILES
