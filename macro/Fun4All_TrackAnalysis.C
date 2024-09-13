//Track + Calo matching for 2024 pp data

#include <G4_ActsGeom.C>
#include <G4_Magnet.C>
#include <GlobalVariables.C>
#include <Trkr_Clustering.C>
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
R__LOAD_LIBRARY(libtrack_reco.so)
R__LOAD_LIBRARY(libtrack_to_calo.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libkfparticle_sphenix.so)

using namespace std;

namespace fs = std::filesystem;

std::string GetFirstLine(std::string listname);
bool is_directory_empty(const fs::path& dir_path);

void Fun4All_TrackAnalysis(
    const int nEvents = 10,
    vector<string> myInputLists = {
        "run46730_0000_trkr_seed.txt",
        "run46730_0000_trkr_cluster.txt",
        "run46730_calo.list"}, 
    bool doTpcOnlyTracking = true,
    bool doEvtDisplay = false,
    bool doEMcalRadiusCorr = true,
    const bool convertSeeds = false)
{
  gSystem->Load("libg4dst.so");

  int verbosity = 0;

  std::cout << "Including " << myInputLists.size() << " files." << std::endl;
  std::string firstfile = GetFirstLine(myInputLists[0]);
  if (*firstfile.c_str() == '\0') return;
  std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment(firstfile);
  int runnumber = runseg.first;
  int segment = runseg.second;

  //The next set of lines figures out folder revisions, file numbers etc
  string outDir = "/sphenix/u/xyu3/hftg01/PhotonConv/macro";

  Enable::DSTOUT = false;
  string outputRecoFileName = "TrackCalo_" + to_string(segment) + "_kfp.root";
  string outputAnaFileName = "TrackCalo_" + to_string(segment) + "_ana.root";
  string outputDstFileName = "TrackCalo_" + to_string(segment) + "_dst.root";

  string outputRecoDir = outDir + "/inReconstruction/" + to_string(runnumber) + "/";
  string makeDirectory = "mkdir -p " + outputRecoDir;
  system(makeDirectory.c_str());
  outputRecoFile = outputRecoDir + outputRecoFileName; // defined in HFReco.h
  string outputAnaFile = outputRecoDir + outputAnaFileName;
  string outputDstFile = outputRecoDir + outputDstFileName;

  string outputJsonDir = outDir + "/inReconstruction/" + to_string(runnumber) + "/EvtDisplay/";
  string makeDirectory2 = "mkdir -p " + outputJsonDir;
  system(makeDirectory2.c_str());

  std::cout << "Reco KFP file: " << outputRecoFile << std::endl;
  std::cout << "Reco ANA file: " << outputAnaFile << std::endl;
  std::cout << "Dst file: " << outputDstFile << std::endl;
  std::cout << "Json file path: " << outputJsonDir << std::endl;

  TpcReadoutInit( runnumber );
  std::cout<< " run: " << runnumber
	   << " samples: " << TRACKING::reco_tpc_maxtime_sample
	   << " pre: " << TRACKING::reco_tpc_time_presample
	   << " vdrift: " << G4TPC::tpc_drift_velocity_reco
	   << std::endl;

  //Create the server
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1);

  //Add all required input files
  for (unsigned int i = 0; i < myInputLists.size(); ++i)
  {
    Fun4AllInputManager *infile = new Fun4AllDstInputManager("DSTin_" + to_string(i));
    std::cout << "Including file " << myInputLists[i] << std::endl;
    //infile->AddFile(myInputLists[i]);
    infile->AddListFile(myInputLists[i]);
    se->registerInputManager(infile);
  }

  auto rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", runnumber);

  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
  rc->set_uint64Flag("TIMESTAMP", runnumber);

  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");

  Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  G4TPC::tpc_drift_velocity_reco = 0.00710;
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
  //to turn on the default static corrections, enable the two lines below
  //G4TPC::ENABLE_STATIC_CORRECTIONS = true;
  //G4TPC::DISTORTIONS_USE_PHI_AS_RADIANS = false;

  G4MAGNET::magfield_rescale = 1;

  TRACKING::pp_mode = false;

  TrackingInit();

  G4TRACKING::convert_seeds_to_svtxtracks = convertSeeds;
  std::cout << "Converting to seeds : " << G4TRACKING::convert_seeds_to_svtxtracks << std::endl;

  /*
   * Either converts seeds to tracks with a straight line/helix fit
   * or run the full Acts track kalman filter fit
   */
  if (G4TRACKING::convert_seeds_to_svtxtracks)
  {
    auto converter = new TrackSeedTrackMapConverter;
    // Option to use TpcTrackSeedContainer or SvtxTrackSeeds
    // can be set to SiliconTrackSeedContainer for silicon-only track fit
    if (doTpcOnlyTracking)
    {
      converter->setTrackSeedName("TpcTrackSeedContainer");
    }
    else
    {
      converter->setTrackSeedName("SvtxTrackSeeds");
    }
    converter->setFieldMap(G4MAGNET::magfield_tracking);
    converter->Verbosity(0);
    se->registerSubsystem(converter);
  }
  else
  {
    auto deltazcorr = new PHTpcDeltaZCorrection;
    deltazcorr->Verbosity(0);
    se->registerSubsystem(deltazcorr);

    // perform final track fit with ACTS
    auto actsFit = new PHActsTrkFitter;
    actsFit->Verbosity(0);
    actsFit->commissioning(G4TRACKING::use_alignment);
    // in calibration mode, fit only Silicons and Micromegas hits
    if (!doTpcOnlyTracking)
    {
      actsFit->fitSiliconMMs(G4TRACKING::SC_CALIBMODE);
      actsFit->setUseMicromegas(G4TRACKING::SC_USE_MICROMEGAS);
    }
    actsFit->set_pp_mode(TRACKING::pp_mode);
    actsFit->set_use_clustermover(true);  // default is true for now
    actsFit->useActsEvaluator(false);
    actsFit->useOutlierFinder(false);
    actsFit->setFieldMap(G4MAGNET::magfield_tracking);
    se->registerSubsystem(actsFit);
  }

  PHSimpleVertexFinder *finder = new PHSimpleVertexFinder;
  finder->Verbosity(0);
  finder->setDcaCut(0.5);
  finder->setTrackPtCut(-99999.);
  finder->setBeamLineCut(1);
  finder->setTrackQualityCut(1000000000);
  if (!doTpcOnlyTracking)
  {
    finder->setRequireMVTX(true);
    finder->setNmvtxRequired(3);
  }
  else
  {
    finder->setRequireMVTX(false);
  }
  finder->setOutlierPairCut(0.1);
  se->registerSubsystem(finder);

  Global_Reco();

  auto projection = new PHActsTrackProjection("CaloProjection");
  float new_cemc_rad = 100.70;//(1-(-0.077))*93.5 recommended cemc radius
  //float new_cemc_rad = 99.1;//(1-(-0.060))*93.5
  //float new_cemc_rad = 97.6;//(1-(-0.044))*93.5, (0.041+0.047)/2=0.044
  if (doEMcalRadiusCorr)
  {
    projection->setLayerRadius(SvtxTrack::CEMC, new_cemc_rad);
  }
  se->registerSubsystem(projection);

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

  //CaloOnly *co = new CaloOnly("CaloOnly", outputRecoFile);
  //se->registerSubsystem(co);

  //TrackOnly *to = new TrackOnly("TracksOnly", outputRecoFile);
  //se->registerSubsystem(to);

  TrackCaloMatch *tcm = new TrackCaloMatch("Tracks_Calo_Match");
  tcm->SetMyTrackMapName("MySvtxTrackMap");
  tcm->writeEventDisplays(doEvtDisplay);
  tcm->setEventDisplayPath(outputJsonDir);
  tcm->setRunDate("2024-06-25");
  tcm->EMcalRadiusUser(doEMcalRadiusCorr);
  tcm->setEMcalRadius(new_cemc_rad);
  tcm->setdphicut(0.2);
  tcm->setdzcut(10000);
  tcm->setTrackPtLowCut(0.5);
  tcm->setEmcalELowCut(0.2);
  se->registerSubsystem(tcm);

  // begin KFParticle
  PhotonConvKFPReco();

  TrackToCalo *ttc = new TrackToCalo("Tracks_And_Calo", outputAnaFile);
  ttc->EMcalRadiusUser(doEMcalRadiusCorr);
  ttc->setEMcalRadius(new_cemc_rad);
  ttc->setKFPtrackMapName("PhotonConv_SvtxTrackMap");
  ttc->setKFPContName("PhotonConv_KFParticle_Container");
  ttc->doTrkrCaloMatching();
  ttc->anaTrkrInfo();
  ttc->anaCaloInfo();
  ttc->setTrackPtLowCut(0.5);
  ttc->setEmcalELowCut(0.2);
  ttc->doTrkrCaloMatching_KFP();
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
    string outputRecoDirMove = outDir + "/Reconstructed/" + to_string(runnumber) + "/";
    string makeDirectoryMove = "mkdir -p " + outputRecoDirMove;
    system(makeDirectoryMove.c_str());
    string moveOutput = "mv " + outputRecoFile + " " + outDir + "/Reconstructed/" + to_string(runnumber);
    std::cout << "moveOutput: " << moveOutput << std::endl;
    system(moveOutput.c_str());
  }

  ifstream file_ana(outputAnaFile.c_str(), ios::binary | ios::ate);
  if (file_ana.good() && (file_ana.tellg() > 100))
  {
    string outputRecoDirMove = outDir + "/Reconstructed/" + to_string(runnumber) + "/";
    string makeDirectoryMove = "mkdir -p " + outputRecoDirMove;
    system(makeDirectoryMove.c_str());
    string moveOutput = "mv " + outputAnaFile + " " + outDir + "/Reconstructed/" + to_string(runnumber);
    std::cout << "moveOutput: " << moveOutput << std::endl;
    system(moveOutput.c_str());
  }

  ifstream file_dst(outputDstFile.c_str(), ios::binary | ios::ate);
  if (file_dst.good() && (file_dst.tellg() > 100))
  {
    string outputRecoDirMove = outDir + "/Reconstructed/" + to_string(runnumber) + "/";
    string makeDirectoryMove = "mkdir -p " + outputRecoDirMove;
    system(makeDirectoryMove.c_str());
    string moveOutput = "mv " + outputDstFile + " " + outDir + "/Reconstructed/" + to_string(runnumber);
    std::cout << "moveOutput: " << moveOutput << std::endl;
    system(moveOutput.c_str());
  }

  fs::path output_json_dir_path = outputJsonDir;
  if (!is_directory_empty(output_json_dir_path))
  {
    string outputJsonDirMove = outDir + "/Reconstructed/" + to_string(runnumber) + "/EvtDisplay/";
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
