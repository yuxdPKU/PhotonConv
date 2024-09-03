input_file=$1

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

  nline_streaming=$(cat $streaminglist | wc -l)
  nline_calo=$(cat $calolist | wc -l)

  echo reading run $runnumber -- nline_streaming: $nline_streaming, nline_calo: $nline_calo
 
done < "$input_file"
