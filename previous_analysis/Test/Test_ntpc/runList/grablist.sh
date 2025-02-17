#!/usr/bin/bash

filename="$1"

#tpc only production 2500 events per DST
#CreateDstList.pl DST_TPC_EVENT_run2pp --cdb 2024p002 --build new --list ${filename}

#streaming production 2500 events per DST (In old runs, 2500 events per DST. In new runs, 100 events per DST. Replacing fast?)
#CreateDstList.pl DST_STREAMING_EVENT_run2pp --cdb 2024p002 --build new --list ${filename}

#streaming fast production 100 events per DST (not used)
#CreateDstList.pl DST_STREAMING_EVENT_run2ppfast --cdb 2024p002 --build new --list ${filename}

#tracking cluster production
#CreateDstList.pl DST_TRKR_CLUSTER_run2pp --cdb 2024p007 --build new --list ${filename}

#tracking seed production
#CreateDstList.pl DST_TRKR_SEED_run2pp --cdb 2024p007 --build new --list ${filename}



#calo production 10000 events per DST in caloy2calib
#CreateDstList.pl DST_CALO_run2pp --cdb 2024p007 --build ana430 --list ${filename}
CreateDstList.pl DST_CALO_run2pp --cdb 2024p007 --build ana435 --list ${filename}

#calo production 10000 events per DST in caloy2test (outdated)
#CreateDstList.pl DST_CALO_run2pp --cdb 2024p006 --build new --list ${filename}
