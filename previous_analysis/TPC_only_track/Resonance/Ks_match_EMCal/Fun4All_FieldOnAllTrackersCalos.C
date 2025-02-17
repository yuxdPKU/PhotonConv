//Track + Calo matching for 2024 pp data
//authors: Xudong Yu <xyu3@bnl.gov>, Antonio Silva <antonio.silva@cern.ch>

#include <G4_ActsGeom.C>
#include <G4_Magnet.C>
#include <GlobalVariables.C>
#include <Trkr_Clustering.C>
#include <Trkr_LaserClustering.C>
#include <Trkr_RecoInit.C>
#include <Trkr_TpcReadoutInit.C>
#include <Trkr_Reco.C>
#include <G4_Global.C>

#include <caloreco/RawClusterBuilderTopo.h>

#include <ffamodules/CDBInterface.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/recoConsts.h>

#include <cdbobjects/CDBTTree.h>

#include <tpccalib/PHTpcResiduals.h>

#include <track_to_calo/TrackCaloMatch.h>
#include <track_to_calo/TrackToCalo.h>
#include <track_to_calo/CaloOnly.h>
#include <track_to_calo/TrackOnly.h>

#include <trackreco/AzimuthalSeeder.h>
#include <trackreco/PHActsTrackProjection.h>
#include <trackingdiagnostics/TrackResiduals.h>

#include <trackbase_historic/SvtxTrack.h>

#include <stdio.h>
#include <iostream>
#include <filesystem>

#include "HFReco.C"

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libcdbobjects.so)
R__LOAD_LIBRARY(libmvtx.so)
R__LOAD_LIBRARY(libintt.so)
R__LOAD_LIBRARY(libtpc.so)
R__LOAD_LIBRARY(libmicromegas.so)
R__LOAD_LIBRARY(libTrackingDiagnostics.so)
R__LOAD_LIBRARY(libtrackingqa.so)
R__LOAD_LIBRARY(libtpcqa.so)
R__LOAD_LIBRARY(libtrack_reco.so)
R__LOAD_LIBRARY(libtrack_to_calo.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libkfparticle_sphenix.so)

using namespace std;

namespace fs = std::filesystem;

std::string GetFirstLine(std::string listname);
bool is_directory_empty(const fs::path& dir_path);

void Fun4All_FieldOnAllTrackersCalos(
    const int nEvents = 10,
    const int runnumber = 46730,
    vector<string> myInputLists = {
        "run46730_0000_trkr.txt",
        "run46730_calo.list"}, 
    std::string outDir = "./",
    bool doTpcOnlyTracking = true,
    bool doEvtDisplay = false,
    bool doEMcalRadiusCorr = true,
    const bool convertSeeds = false)
{
  gSystem->Load("libg4dst.so");

  int verbosity = 0;

  G4TRACKING::convert_seeds_to_svtxtracks = convertSeeds;
  std::cout << "Converting to seeds : " << G4TRACKING::convert_seeds_to_svtxtracks << std::endl;
  std::cout << "Including " << myInputLists.size() << " files." << std::endl;
  std::string firstfile = GetFirstLine(myInputLists[0]);
  if (*firstfile.c_str() == '\0') return;
  std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment(firstfile);
  //int runnumber = runseg.first;
  //int segment = runseg.second;
  int segment = 0;

 TRACKING::pp_mode = true;

  // distortion calibration mode
  /*
   * set to true to enable residuals in the TPC with
   * TPC clusters not participating to the ACTS track fit
   */
  G4TRACKING::SC_CALIBMODE = false;

  ACTSGEOM::mvtxMisalignment = 100;
  ACTSGEOM::inttMisalignment = 100.;
  ACTSGEOM::tpotMisalignment = 100.;

  Enable::DSTOUT = false;

  string outputRecoFileName = "TrackCalo_run" + to_string(runnumber) + "_kfp.root";
  string outputValidateFileName = "TrackCalo_run" + to_string(runnumber) + "_ks.root";
  string outputAnaFileName = "TrackCalo_run" + to_string(runnumber) + "_ana.root";
  string outputDstFileName = "TrackCalo_run" + to_string(runnumber) + "_dst.root";

  string outputRecoDir = outDir + "/inReconstruction/";
  string makeDirectory = "mkdir -p " + outputRecoDir;
  system(makeDirectory.c_str());
  outputRecoFile = outputRecoDir + outputRecoFileName; // defined in HFReco.h
  outputValidateFile = outputRecoDir + outputValidateFileName; // defined in HFReco.h
  string outputAnaFile = outputRecoDir + outputAnaFileName;
  string outputDstFile = outputRecoDir + outputDstFileName;

  string outputJsonDir = outDir + "/inReconstruction/EvtDisplay/";
  string makeDirectory2 = "mkdir -p " + outputJsonDir;
  system(makeDirectory2.c_str());

  std::cout << "Reco KFP file: " << outputRecoFile << std::endl;
  std::cout << "Reco Ks Validation file: " << outputValidateFile << std::endl;
  std::cout << "Reco ANA file: " << outputAnaFile << std::endl;
  std::cout << "Dst file: " << outputDstFile << std::endl;
  std::cout << "Json file path: " << outputJsonDir << std::endl;

  //Create the server
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(2);
  auto rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", runnumber);
  rc->set_IntFlag("RUNSEGMENT", segment);

  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
  rc->set_uint64Flag("TIMESTAMP", runnumber);
  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");

  Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  TpcReadoutInit( runnumber );
  std::cout<< " run: " << runnumber
	   << " samples: " << TRACKING::reco_tpc_maxtime_sample
	   << " pre: " << TRACKING::reco_tpc_time_presample
	   << " vdrift: " << G4TPC::tpc_drift_velocity_reco
	   << std::endl;

  CDBInterface *cdb = CDBInterface::instance();
  std::string tpc_dv_calib_dir = cdb->getUrl("TPC_DRIFT_VELOCITY");
  if (tpc_dv_calib_dir.empty())
  {
    std::cout << "No calibrated TPC drift velocity for Run " << runnumber << ". Use default value " << G4TPC::tpc_drift_velocity_reco << " cm/ns" << std::endl;
  }
  else
  {
    CDBTTree *cdbttree = new CDBTTree(tpc_dv_calib_dir);
    cdbttree->LoadCalibrations();
    G4TPC::tpc_drift_velocity_reco = cdbttree->GetSingleFloatValue("tpc_drift_velocity");
    std::cout << "Use calibrated TPC drift velocity for Run " << runnumber << ": " << G4TPC::tpc_drift_velocity_reco << " cm/ns" << std::endl;
  }

  G4TPC::ENABLE_MODULE_EDGE_CORRECTIONS = true;
  // enable TPC zero suppression
  if(runnumber>51428)
  {
    TRACKING::tpc_zero_supp = true;
  }

  //to turn on the default static corrections, enable the two lines below
  //G4TPC::ENABLE_STATIC_CORRECTIONS = true;
  //G4TPC::USE_PHI_AS_RAD_STATIC_CORRECTIONS = false;

  //to turn on the average corrections derived from simulation, enable the three lines below
  //note: these are designed to be used only if static corrections are also applied
  //G4TPC::ENABLE_AVERAGE_CORRECTIONS = true;
  //G4TPC::USE_PHI_AS_RAD_AVERAGE_CORRECTIONS = false;
  //G4TPC:average_correction_filename = std::string(getenv("CALIBRATIONROOT")) + "/distortion_maps/average_minus_static_distortion_inverted_10-new.root";

  G4MAGNET::magfield_rescale = 1;
  TrackingInit();

  //Add all required input files
  for (unsigned int i = 0; i < myInputLists.size(); ++i)
  {
    Fun4AllInputManager *infile = new Fun4AllDstInputManager("DSTin_" + to_string(i));
    std::cout << "Including file " << myInputLists[i] << std::endl;
    //infile->AddFile(myInputLists[i]);
    infile->AddListFile(myInputLists[i]);
    se->registerInputManager(infile);
  }

  Global_Reco();

  auto projection = new PHActsTrackProjection("CaloProjection");
  float new_cemc_rad = 100.70;//(1-(-0.077))*93.5 recommended cemc radius at shower max
  //float new_cemc_rad = 99.1;//(1-(-0.060))*93.5
  //float new_cemc_rad = 97.6;//(1-(-0.044))*93.5, (0.041+0.047)/2=0.044
  if (doEMcalRadiusCorr)
  {
    projection->setLayerRadius(SvtxTrack::CEMC, new_cemc_rad);
  }
  se->registerSubsystem(projection);

  std::cout << "Building clusters" << std::endl;
  RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
  ClusterBuilder->Detector("CEMC");
  ClusterBuilder->set_threshold_energy(0.030);  // for when using basic calibration
  std::string emc_prof = getenv("CALIBRATIONROOT");
  emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
  ClusterBuilder->LoadProfile(emc_prof);
  ClusterBuilder->set_UseTowerInfo(1);  // to use towerinfo objects rather than old RawTower
  se->registerSubsystem(ClusterBuilder);

/*
  RawClusterBuilderTopo* ClusterBuilder1 = new RawClusterBuilderTopo("EMcalRawClusterBuilderTopo1");
  ClusterBuilder1->Verbosity(verbosity);
  ClusterBuilder1->set_nodename("TOPOCLUSTER_EMCAL");
  ClusterBuilder1->set_enable_HCal(false);
  ClusterBuilder1->set_enable_EMCal(true);
  //ClusterBuilder1->set_noise(0.0025, 0.006, 0.03);
  ClusterBuilder1->set_noise(0.01, 0.03, 0.03);
  ClusterBuilder1->set_significance(4.0, 2.0, 1.0);
  ClusterBuilder1->allow_corner_neighbor(true);
  ClusterBuilder1->set_do_split(true);
  ClusterBuilder1->set_minE_local_max(1.0, 2.0, 0.5);
  ClusterBuilder1->set_R_shower(0.025);
  se->registerSubsystem(ClusterBuilder1);
*/

  //For particle flow studies
  RawClusterBuilderTopo* ClusterBuilder2 = new RawClusterBuilderTopo("EMcalRawClusterBuilderTopo2");
  ClusterBuilder2->Verbosity(verbosity);
  ClusterBuilder2->set_nodename("TOPOCLUSTER_HCAL");
  ClusterBuilder2->set_enable_HCal(true);
  ClusterBuilder2->set_enable_EMCal(false);
  //ClusterBuilder2->set_noise(0.0025, 0.006, 0.03);
  ClusterBuilder2->set_noise(0.01, 0.03, 0.03);
  ClusterBuilder2->set_significance(4.0, 2.0, 1.0);
  ClusterBuilder2->allow_corner_neighbor(true);
  ClusterBuilder2->set_do_split(true);
  ClusterBuilder2->set_minE_local_max(1.0, 2.0, 0.5);
  ClusterBuilder2->set_R_shower(0.025);
  se->registerSubsystem(ClusterBuilder2);

  TrackCaloMatch *tcm = new TrackCaloMatch("Tracks_Calo_Match");
  tcm->SetMyTrackMapName("MySvtxTrackMap");
  tcm->writeEventDisplays(doEvtDisplay);
  tcm->setEventDisplayPath(outputJsonDir);
  tcm->setRunDate("2024-06-25");
  tcm->EMcalRadiusUser(doEMcalRadiusCorr);
  tcm->setEMcalRadius(new_cemc_rad);
  tcm->setdphicut(0.15);
  tcm->setdzcut(10);
  tcm->setTrackPtLowCut(0.2);
  tcm->setEmcalELowCut(0.0);
  tcm->setnTpcClusters(20);
  tcm->setTrackQuality(200);
  tcm->setRawClusContEMName("CLUSTERINFO_CEMC");
  se->registerSubsystem(tcm);

  // begin KFParticle
  KsReco(); // Ks reco for validation

  TrackToCalo *ttc = new TrackToCalo("Tracks_And_Calo", outputAnaFile);
  ttc->EMcalRadiusUser(doEMcalRadiusCorr);
  ttc->setEMcalRadius(new_cemc_rad);
  ttc->setKFPtrackMapName("PhotonConv_SvtxTrackMap");
  ttc->setKFPContName("PhotonConv_KFParticle_Container");
  //ttc->doTrkrCaloMatching();
  //ttc->anaTrkrInfo();
  //ttc->anaCaloInfo();
  ttc->setTrackPtLowCut(0.2);
  ttc->setEmcalELowCut(0.0);
  ttc->setnTpcClusters(20);
  ttc->setTrackQuality(200);
  ttc->doTrkrCaloMatching_KFP();
  ttc->setRawClusContEMName("CLUSTERINFO_CEMC");
  se->registerSubsystem(ttc);

  if (Enable::DSTOUT)
  {
    Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", outputDstFile);
    //out->StripNode("RUN");
    //out->AddNode("Sync");
    out->SaveRunNode(0);
    se->registerOutputManager(out);
  }

  se->run(nEvents);
  se->End();
  se->PrintTimer();

  ifstream file_reco(outputRecoFile.c_str(), ios::binary | ios::ate);
  if (file_reco.good() && (file_reco.tellg() > 100))
  {
    string outputRecoDirMove = outDir + "/Reconstructed/";
    string makeDirectoryMove = "mkdir -p " + outputRecoDirMove;
    system(makeDirectoryMove.c_str());
    string moveOutput = "mv " + outputRecoFile + " " + outDir + "/Reconstructed/";
    std::cout << "moveOutput: " << moveOutput << std::endl;
    system(moveOutput.c_str());
  }

  ifstream file_validation(outputValidateFile.c_str(), ios::binary | ios::ate);
  if (file_validation.good() && (file_validation.tellg() > 100))
  {
    string outputRecoDirMove = outDir + "/Reconstructed/";
    string makeDirectoryMove = "mkdir -p " + outputRecoDirMove;
    system(makeDirectoryMove.c_str());
    string moveOutput = "mv " + outputValidateFile + " " + outDir + "/Reconstructed/";
    std::cout << "moveOutput: " << moveOutput << std::endl;
    system(moveOutput.c_str());
  }

  ifstream file_ana(outputAnaFile.c_str(), ios::binary | ios::ate);
  if (file_ana.good() && (file_ana.tellg() > 100))
  {
    string outputRecoDirMove = outDir + "/Reconstructed/";
    string makeDirectoryMove = "mkdir -p " + outputRecoDirMove;
    system(makeDirectoryMove.c_str());
    string moveOutput = "mv " + outputAnaFile + " " + outDir + "/Reconstructed/";
    std::cout << "moveOutput: " << moveOutput << std::endl;
    system(moveOutput.c_str());
  }

  ifstream file_dst(outputDstFile.c_str(), ios::binary | ios::ate);
  if (file_dst.good() && (file_dst.tellg() > 100))
  {
    string outputRecoDirMove = outDir + "/Reconstructed/";
    string makeDirectoryMove = "mkdir -p " + outputRecoDirMove;
    system(makeDirectoryMove.c_str());
    string moveOutput = "mv " + outputDstFile + " " + outDir + "/Reconstructed/";
    std::cout << "moveOutput: " << moveOutput << std::endl;
    system(moveOutput.c_str());
  }

  fs::path output_json_dir_path = outputJsonDir;
  if (!is_directory_empty(output_json_dir_path))
  {
    string outputJsonDirMove = outDir + "/Reconstructed/EvtDisplay/";
    string makeDirectoryMove = "mkdir -p " + outputJsonDirMove;
    system(makeDirectoryMove.c_str());
    string moveOutput = "mv " + outputJsonDir + "* " + outputJsonDirMove;
    std::cout << "moveOutput: " << moveOutput << std::endl;
    system(moveOutput.c_str());
  }

  delete se;
  std::cout << "All done" << std::endl;
  gSystem->Exit(0);

  return;
}

std::string GetFirstLine(std::string listname)
{
  std::ifstream file(listname);

  std::string firstLine = "";
  if (file.is_open()) {
      if (std::getline(file, firstLine)) {
          std::cout << "First Line: " << firstLine << std::endl;
      } else {
          std::cerr << "Unable to read first line of file" << std::endl;
      }
      file.close();
  } else {
      std::cerr << "Unable to open file" << std::endl;
  }
  return firstLine;
}

bool is_directory_empty(const fs::path& dir_path) {
    if (fs::exists(dir_path) && fs::is_directory(dir_path)) {
        return fs::directory_iterator(dir_path) == fs::directory_iterator();
    }
    return false;
}
