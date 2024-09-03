input_file=$1

if [ ! -f "$input_file" ]; then
  echo "File not found: $input_file"
  exit 1
fi


while IFS= read -r line; do
  read runnumber jobnumber <<< "$line"
  echo Run $runnumber, Njobs $jobnumber
  cp condor-data-trkrcalo-template.job condor-data-trkrcalo-$runnumber.job
  sed -i "s/RUNNUMBER/$runnumber/g" condor-data-trkrcalo-$runnumber.job
  sed -i "s/JOBNUMBER/$jobnumber/g" condor-data-trkrcalo-$runnumber.job
  condor_submit condor-data-trkrcalo-$runnumber.job
done < "$input_file"
