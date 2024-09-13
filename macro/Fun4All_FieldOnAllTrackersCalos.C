//Track + Calo matching for 2024 pp data
//authors: Xudong Yu <xyu3@bnl.gov>, Antonio Silva <antonio.silva@cern.ch>

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

void Fun4All_FieldOnAllTrackersCalos(
    const int nEvents = 10,
    vector<string> myInputLists = {
        "run46730_0000_trkr.txt",
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

  G4TRACKING::convert_seeds_to_svtxtracks = convertSeeds;
  std::cout << "Converting to seeds : " << G4TRACKING::convert_seeds_to_svtxtracks << std::endl;

  // enable TPC zero suppression
  if(runnumber>51428)
  {
    TRACKING::tpc_zero_supp = true;
  }

  TrackingInit();
  if(doTpcOnlyTracking)
  {
    Tpc_HitUnpacking();
  }
  else
  {
    Mvtx_HitUnpacking();
    Intt_HitUnpacking();
    Tpc_HitUnpacking();
    Micromegas_HitUnpacking();
  }

  Mvtx_Clustering();
  Intt_Clustering();

  auto tpcclusterizer = new TpcClusterizer;
  tpcclusterizer->Verbosity(verbosity);
  tpcclusterizer->set_do_hit_association(G4TPC::DO_HIT_ASSOCIATION);
  tpcclusterizer->set_rawdata_reco();
  se->registerSubsystem(tpcclusterizer);

  Micromegas_Clustering();

  /*
   * Begin Track Seeding
   */

  /*
   * Silicon Seeding
   */
  auto silicon_Seeding = new PHActsSiliconSeeding;
  silicon_Seeding->Verbosity(0);
  silicon_Seeding->searchInIntt();
  silicon_Seeding->setinttRPhiSearchWindow(0.4);
  silicon_Seeding->setinttZSearchWindow(1.6);
  silicon_Seeding->seedAnalysis(false);
  se->registerSubsystem(silicon_Seeding);

  auto merger = new PHSiliconSeedMerger;
  merger->Verbosity(0);
  se->registerSubsystem(merger);

  /*
   * Tpc Seeding
   */
  auto seeder = new PHCASeeding("PHCASeeding");
  double fieldstrength = std::numeric_limits<double>::quiet_NaN();  // set by isConstantField if constant
  bool ConstField = isConstantField(G4MAGNET::magfield_tracking, fieldstrength);
  if (ConstField)
  {
    seeder->useConstBField(true);
    seeder->constBField(fieldstrength);
  }
  else
  {
    seeder->set_field_dir(-1 * G4MAGNET::magfield_rescale);
    seeder->useConstBField(false);
    seeder->magFieldFile(G4MAGNET::magfield_tracking);  // to get charge sign right
  }
  seeder->Verbosity(verbosity);
  seeder->SetLayerRange(7, 55);
//used in https://github.com/klendathu2k/slurp-examples/blob/87e186c1149cccd0002ee2cd2f5cad8ceaa362e9/sPHENIX/TrackingProduction/Fun4All_JobA.C#L120
  seeder->SetSearchWindow(1.5,0.05); // z-width and phi-width
//used in https://github.com/sPHENIX-Collaboration/macros/blob/b9d924261ee96352f26867bc4bf80fd37d0a2d04/TrackingProduction/Fun4All_FieldOnAllTrackers.C#L191
//  seeder->SetSearchWindow(2.,0.05); // z-width and phi-width, default in macro at 1.5 and 0.05
//  seeder->SetClusAdd_delta_window(3.0,0.06); //  (0.5, 0.005) are default; sdzdr_cutoff, d2/dr2(phi)_cutoff
//  //seeder->SetNClustersPerSeedRange(4,60); // default is 6, 6
  seeder->SetMinHitsPerCluster(0);
  seeder->SetMinClustersPerTrack(3);
  seeder->useFixedClusterError(true);
  seeder->set_pp_mode(TRACKING::pp_mode);
  se->registerSubsystem(seeder);

  // expand stubs in the TPC using simple kalman filter
  auto cprop = new PHSimpleKFProp("PHSimpleKFProp");
  cprop->set_field_dir(G4MAGNET::magfield_rescale);
  if (ConstField)
  {
    cprop->useConstBField(true);
    cprop->setConstBField(fieldstrength);
  }
  else
  {
    cprop->magFieldFile(G4MAGNET::magfield_tracking);
    cprop->set_field_dir(-1 * G4MAGNET::magfield_rescale);
  }
  cprop->useFixedClusterError(true);
  cprop->set_max_window(5.);
  cprop->Verbosity(verbosity);
  cprop->set_pp_mode(TRACKING::pp_mode);
  se->registerSubsystem(cprop);

  /*
   * Track Matching between silicon and TPC
   */
  // The normal silicon association methods
  // Match the TPC track stubs from the CA seeder to silicon track stubs from PHSiliconTruthTrackSeeding
  auto silicon_match = new PHSiliconTpcTrackMatching;
  silicon_match->Verbosity(0);
  silicon_match->set_x_search_window(2.);
  silicon_match->set_y_search_window(2.);
  silicon_match->set_z_search_window(5.);
  silicon_match->set_phi_search_window(0.2);
  silicon_match->set_eta_search_window(0.1);
  silicon_match->set_use_old_matching(true);
  silicon_match->set_pp_mode(true);
  se->registerSubsystem(silicon_match);

  // Match TPC track stubs from CA seeder to clusters in the micromegas layers
  auto mm_match = new PHMicromegasTpcTrackMatching;
  mm_match->Verbosity(0);
  mm_match->set_rphi_search_window_lyr1(0.4);
  mm_match->set_rphi_search_window_lyr2(13.0);
  mm_match->set_z_search_window_lyr1(26.0);
  mm_match->set_z_search_window_lyr2(0.4);

  mm_match->set_min_tpc_layer(38);             // layer in TPC to start projection fit
  mm_match->set_test_windows_printout(false);  // used for tuning search windows only
  se->registerSubsystem(mm_match);

  /*
   * End Track Seeding
   */

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
