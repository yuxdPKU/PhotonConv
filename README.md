# PhotonConv
Photon conversion dielectron analysis for sPHENIX 2024 pp data

- Compile analysis module
```
cd module
mkdir build
cd build
../autogen.sh --prefix=$MYINSTALL
make -j 4
make install
```

- Make DST list, for example 6x6 bunch crossing data
```
cd macro
./make_tpclist_6x6.sh
```

- `run[runnumber]_tpc.list` and `run[runnumber]_calo.list` will show up in runList/ directory

- Each TPC DST has 2500 events (fast type of DST only has 100 events), split each DST into one file in trackrunlist folder
```
./splitlist.sh [runnumber]
```

- Submit jobs in condorJob
**Do not forget to change Initialdir, RunNumber, input DST list in Arguments, DVCorrTag which is related with pre-calib drift velocity, and the number of jobs you want to submit**
```
Initialdir     = /sphenix/u/xyu3/hftg01
Executable     = $(Initialdir)/run_data.sh
RunNumber      = 51103
DVCorrTag      = 2
Arguments      = "./runList/trackrunlist/run$(RunNumber)_$INT(Process,%04d).txt ./runList/run$(RunNumber)_calo.list $(DVCorrTag)"
Queue 9
```

- Output root file saved in Reconstructed/[runnumber]

- Offline analysis code in OfflineAna
