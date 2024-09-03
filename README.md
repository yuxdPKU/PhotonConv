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

- Prepare a run list, one line one runnumber, e.g. `macro/runList/run51619-51886_1000seg_condor_trkrcalo.list`

- Grab DST file list by using CreateDstList.pl:
```
cd macro/runList
./grablist.sh run51619-51886_1000seg_condor_trkrcalo.list
```
different kinds of DST list are generated: streaming production with all trackers, streaming production with tpc only, tracking production in cluster level, tracking production in seed level, calo production

- Split tracker DST list into several files
```
./splitlist_list_trkr_calo.sh [in-file-list-name] [out-file-list-name] [number-of-DST-per-file]
```
or
```
./splitlist_list_tpc.sh [in-file-list-name] [out-file-list-name] [number-of-DST-per-file]
```

- Submit jobs in condorJob
If you want to do tracking only with track fitting (tracking production as input), use `gen_condor_job_trkrcalo.sh`.
```
./gen_condor_job.sh [path-of-list-file]
```
If you want to do tracking from hit unpacker (streaming production as input), use `gen_condor_job.sh`.
```
./gen_condor_job_trkrcalo.sh [path-of-list-file]
```
**Do not forget to change Initialdir!!!**

- When job running, output root file temporarily stored in inReconstruction/[runnumber]

- When job finished, output root file saved in Reconstructed/[runnumber]

- Offline analysis code in OfflineAna
