/*
 * This macro shows a minimum working example of running the tracking
 * hit unpackers with some basic seeding algorithms to try to put together
 * tracks. There are some analysis modules run at the end which package
 * hits, clusters, and clusters on tracks into trees for analysis.
 */

#include <fun4all/Fun4AllUtils.h>
#include <G4_ActsGeom.C>
#include <G4_Global.C>
#include <G4_Magnet.C>
#include <G4_Mbd.C>
#include <GlobalVariables.C>
#include <QA.C>
#include <Calo_Calib.C>
#include <Trkr_Clustering.C>
#include <Trkr_LaserClustering.C>
#include <Trkr_Reco.C>
#include <Trkr_RecoInit.C>
#include <Trkr_TpcReadoutInit.C>

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/recoConsts.h>

#include <cdbobjects/CDBTTree.h>

#include <tpccalib/PHTpcResiduals.h>

#include <trackingqa/InttClusterQA.h>
#include <trackingqa/MicromegasClusterQA.h>
#include <trackingqa/MvtxClusterQA.h>
#include <trackingqa/TpcClusterQA.h>
#include <tpcqa/TpcRawHitQA.h>
#include <trackingqa/TpcSeedsQA.h>

#include <trackreco/PHActsTrackProjection.h>

#include <trackingdiagnostics/TrackResiduals.h>
#include <trackingdiagnostics/TrkrNtuplizer.h>
//#include <trackingdiagnostics/KshortReconstruction.h>

#include <track_to_calo/TrackCaloMatch.h>
#include <track_to_calo/TrackToCalo.h>
#include <track_to_calo/CaloOnly.h>
#include <track_to_calo/TrackOnly.h>

#include <caloreco/CaloGeomMapping.h>
#include <caloreco/CaloGeomMappingv2.h>
#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/RawClusterBuilderTopo.h>

#include <stdio.h>

#pragma GCC diagnostic push

#pragma GCC diagnostic ignored "-Wundefined-internal"

#include <kfparticle_sphenix/KFParticle_sPHENIX.h>
#include <kfparticle_sphenix/KshortReconstruction_local.h>

#pragma GCC diagnostic pop

void KFPReco(std::string module_name = "KFPReco", std::string decaydescriptor = "K_S0 -> pi^+ pi^-", std::string outfile = "KFP.root", std::string trackmapName = "SvtxTrackMap", std::string containerName = "KFParticle");

R__LOAD_LIBRARY(libkfparticle_sphenix.so)
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
R__LOAD_LIBRARY(libtrack_to_calo.so)
R__LOAD_LIBRARY(libtrack_reco.so)
R__LOAD_LIBRARY(libcalo_reco.so)
//R__LOAD_LIBRARY(libcalogeomtest.so)
R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(libcentrality.so)
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libepd.so)
R__LOAD_LIBRARY(libzdcinfo.so)
void Fun4All_TrackFitting_PhotonConv(
    const int nEvents = 10,
    const std::string clusterfilename = "DST_TRKR_CLUSTER_run2pp_ana466_2024p012_v001-00053534-00000.root",
    const std::string clusterdir = "/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana466_2024p012_v001/DST_TRKR_CLUSTER/run_00053500_00053600/dst/",
    const std::string seedfilename = "DST_TRKR_SEED_run2pp_ana468_2024p012_v002-00053534-00000.root",
    const std::string seeddir = "/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana468_2024p012_v002/DST_TRKR_SEED/run_00053500_00053600/dst/",
    const std::string calofilename = "DST_CALO_run2pp_ana462_2024p010_v001-00053534-00000.root",
    const std::string calodir = "/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana462_2024p010_v001/DST_CALO/run_00053500_00053600/dst/",
    const std::string outfilename = "clusters_seeds",
    const std::string outdir = "./root",
    const int index = 0,
    const int stepsize = 10,
    const bool convertSeeds = false)
{
  std::string inputTrkrClusterFile = clusterdir + clusterfilename;
  std::string inputTrkrSeedFile = seeddir + seedfilename;
  std::string inputCaloFile = calodir + calofilename;

  G4TRACKING::convert_seeds_to_svtxtracks = convertSeeds;
  std::cout << "Converting to seeds : " << G4TRACKING::convert_seeds_to_svtxtracks << std::endl;
  std::pair<int, int>
      runseg = Fun4AllUtils::GetRunSegment(seedfilename);
  int runnumber = runseg.first;
  int segment = runseg.second;

  string outDir = outdir + "/inReconstruction/" + to_string(runnumber) + "/";
  string makeDirectory = "mkdir -p " + outDir;
  system(makeDirectory.c_str());
  TString outfile = outDir + outfilename + "_" + runnumber + "-" + segment + "-" + index + ".root";
  std::cout<<"outfile "<<outfile<<std::endl;
  std::string theOutfile = outfile.Data();

  G4TRACKING::SC_CALIBMODE = false;
  Enable::MVTX_APPLYMISALIGNMENT = true;
  ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT;
  TRACKING::pp_mode = true;

  auto se = Fun4AllServer::instance();
  se->Verbosity(1);
  auto rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", runnumber);
  rc->set_IntFlag("RUNSEGMENT", segment);

  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");
  rc->set_uint64Flag("TIMESTAMP", runnumber);
  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");

  TpcReadoutInit(runnumber);
 // these lines show how to override the drift velocity and time offset values set in TpcReadoutInit
  // G4TPC::tpc_drift_velocity_reco = 0.0073844; // cm/ns
  // TpcClusterZCrossingCorrection::_vdrift = G4TPC::tpc_drift_velocity_reco;
  // G4TPC::tpc_tzero_reco = -5*50;  // ns
  std::cout << " run: " << runnumber
            << " samples: " << TRACKING::reco_tpc_maxtime_sample
            << " pre: " << TRACKING::reco_tpc_time_presample
            << " vdrift: " << G4TPC::tpc_drift_velocity_reco
            << std::endl;

  // distortion calibration mode
  /*
   * set to true to enable residuals in the TPC with
   * TPC clusters not participating to the ACTS track fit
   */
  G4TRACKING::SC_CALIBMODE = false;
  TRACKING::pp_mode = true;

  ACTSGEOM::mvtxMisalignment = 100;
  ACTSGEOM::inttMisalignment = 100.;
  ACTSGEOM::tpotMisalignment = 100.;


  Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  G4TPC::ENABLE_MODULE_EDGE_CORRECTIONS = true;

  //to turn on the default static corrections, enable the two lines below
  G4TPC::ENABLE_STATIC_CORRECTIONS = true;
  G4TPC::USE_PHI_AS_RAD_STATIC_CORRECTIONS = false;

  //to turn on the average corrections derived from simulation, enable the three lines below
  //note: these are designed to be used only if static corrections are also applied
  G4TPC::ENABLE_AVERAGE_CORRECTIONS = true;
  G4TPC::USE_PHI_AS_RAD_AVERAGE_CORRECTIONS = false;
  G4TPC::average_correction_filename = "/sphenix/tg/tg01/jets/bkimelman/BenProduction/Feb25_2025/Laminations_run2pp_ana466_2024p012_v001-000" + to_string(runnumber) + ".root";
  std::cout<<"Average distortion map used: "<<G4TPC::average_correction_filename<<std::endl;
  //G4TPC::average_correction_filename = std::string(getenv("CALIBRATIONROOT")) + "/distortion_maps/average_minus_static_distortion_inverted_10-new.root";

  G4MAGNET::magfield_rescale = 1;
  TrackingInit();

  auto hitsin_seed = new Fun4AllDstInputManager("DSTin_seed");
  hitsin_seed->fileopen(inputTrkrSeedFile);
  se->registerInputManager(hitsin_seed);

  auto hitsin_cluster = new Fun4AllDstInputManager("DSTin_cluster");
  hitsin_cluster->fileopen(inputTrkrClusterFile);
  se->registerInputManager(hitsin_cluster);

  auto hitsin_calo = new Fun4AllDstInputManager("DSTin_calo");
  hitsin_calo->fileopen(inputCaloFile);
  se->registerInputManager(hitsin_calo);

  G4TPC::REJECT_LASER_EVENTS = true;
  // reject laser events if G4TPC::REJECT_LASER_EVENTS is true
  Reject_Laser_Events();

  // Always apply preliminary distortion corrections to TPC clusters before silicon matching
  // and refit the trackseeds. Replace KFProp fits with the new fit parameters in the TPC seeds.
  auto prelim_distcorr = new PrelimDistortionCorrection;
  prelim_distcorr->set_pp_mode(true);
  prelim_distcorr->Verbosity(0);
  se->registerSubsystem(prelim_distcorr);

  /*
   * Track Matching between silicon and TPC
   */
  // The normal silicon association methods
  // Match the TPC track stubs from the CA seeder to silicon track stubs from PHSiliconTruthTrackSeeding
  auto silicon_match = new PHSiliconTpcTrackMatching;
  silicon_match->Verbosity(0);
  silicon_match->set_use_legacy_windowing(false);
  silicon_match->set_pp_mode(TRACKING::pp_mode);
  if(G4TPC::ENABLE_AVERAGE_CORRECTIONS)
    {
      // reset phi matching window to be centered on zero
      // it defaults to being centered on -0.1 radians for the case of static corrections only
      std::array<double,3> arrlo = {-0.15,0,0};
      std::array<double,3> arrhi = {0.15,0,0};
      silicon_match->window_dphi.set_QoverpT_range(arrlo, arrhi);
    }
    se->registerSubsystem(silicon_match);

  // Match TPC track stubs from CA seeder to clusters in the micromegas layers
  auto mm_match = new PHMicromegasTpcTrackMatching;
  mm_match->Verbosity(0);
  mm_match->set_pp_mode(TRACKING::pp_mode);
  mm_match->set_rphi_search_window_lyr1(3.);
  mm_match->set_rphi_search_window_lyr2(15.0);
  mm_match->set_z_search_window_lyr1(30.0);
  mm_match->set_z_search_window_lyr2(3.);

  mm_match->set_min_tpc_layer(38);             // layer in TPC to start projection fit
  mm_match->set_test_windows_printout(false);  // used for tuning search windows only
  se->registerSubsystem(mm_match);


  /*
   * Either converts seeds to tracks with a straight line/helix fit
   * or run the full Acts track kalman filter fit
   */
  if (G4TRACKING::convert_seeds_to_svtxtracks)
  {
    auto converter = new TrackSeedTrackMapConverter;
    // Default set to full SvtxTrackSeeds. Can be set to
    // SiliconTrackSeedContainer or TpcTrackSeedContainer
    converter->setTrackSeedName("SvtxTrackSeedContainer");
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
    actsFit->fitSiliconMMs(G4TRACKING::SC_CALIBMODE);
    actsFit->setUseMicromegas(G4TRACKING::SC_USE_MICROMEGAS);
    actsFit->set_pp_mode(TRACKING::pp_mode);
    actsFit->set_use_clustermover(true);  // default is true for now
    actsFit->useActsEvaluator(false);
    actsFit->useOutlierFinder(false);
    actsFit->setFieldMap(G4MAGNET::magfield_tracking);
    se->registerSubsystem(actsFit);

    auto cleaner = new PHTrackCleaner();
    cleaner->Verbosity(0);
    cleaner->set_pp_mode(TRACKING::pp_mode);
    se->registerSubsystem(cleaner);

    if (G4TRACKING::SC_CALIBMODE)
    {
      /*
       * in calibration mode, calculate residuals between TPC and fitted tracks,
       * store in dedicated structure for distortion correction
       */
      auto residuals = new PHTpcResiduals;
      const TString tpc_residoutfile = theOutfile + "_PhTpcResiduals.root";
      residuals->setOutputfile(tpc_residoutfile.Data());
      residuals->setUseMicromegas(G4TRACKING::SC_USE_MICROMEGAS);

      // matches Tony's analysis
      residuals->setMinPt(0.2);

      // reconstructed distortion grid size (phi, r, z)
      residuals->setGridDimensions(36, 48, 80);
      se->registerSubsystem(residuals);
    }
  }

  auto finder = new PHSimpleVertexFinder;
  finder->Verbosity(0);
  finder->setDcaCut(0.5);
  finder->setTrackPtCut(-99999.);
  finder->setBeamLineCut(1);
  finder->setTrackQualityCut(1000000000);
  finder->setNmvtxRequired(3);
  finder->setOutlierPairCut(0.1);
  se->registerSubsystem(finder);

  Global_Reco();

  bool doEMcalRadiusCorr = true;
  auto projection = new PHActsTrackProjection("CaloProjection");
  float new_cemc_rad = 99.; // from DetailedCalorimeterGeometry, project to inner surface
  //float new_cemc_rad = 100.70;//(1-(-0.077))*93.5 recommended cemc radius at shower max
  //float new_cemc_rad = 99.1;//(1-(-0.060))*93.5
  //float new_cemc_rad = 97.6;//(1-(-0.044))*93.5, (0.041+0.047)/2=0.044
  if (doEMcalRadiusCorr)
  {
    projection->setLayerRadius(SvtxTrack::CEMC, new_cemc_rad);
  }
  se->registerSubsystem(projection);

  /////////////////////////////////////////////////////
  // Set status of CALO towers, Calibrate towers,  Cluster
  //Process_Calo_Calib();

  //my calo reco
  std::cout<<"Begin my calo reco"<<std::endl;
  // Load the modified geometry
  CaloGeomMappingv2 *cgm = new CaloGeomMappingv2();
  cgm->set_detector_name("CEMC");
  cgm->setTowerGeomNodeName("TOWERGEOM_CEMCv3");
  se->registerSubsystem(cgm);

  //////////////////
  // Clusters
  std::cout << "Building clusters" << std::endl;
  RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
  ClusterBuilder->Detector("CEMC");
  ClusterBuilder->setUseRawTowerGeomv5(true);
  ClusterBuilder->setProjectToInnerSurface(true);
  ClusterBuilder->set_threshold_energy(0.070);  // for when using basic calibration
  std::string emc_prof = getenv("CALIBRATIONROOT");
  emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
  ClusterBuilder->LoadProfile(emc_prof);
  ClusterBuilder->set_UseTowerInfo(1);  // to use towerinfo objects rather than old RawTower
  se->registerSubsystem(ClusterBuilder);

  //For particle flow studies
  RawClusterBuilderTopo* ClusterBuilder2 = new RawClusterBuilderTopo("EMcalRawClusterBuilderTopo2");
  ClusterBuilder2->Verbosity(0);
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
  tcm->writeEventDisplays(false);
  tcm->EMcalRadiusUser(doEMcalRadiusCorr);
  tcm->setEMcalRadius(new_cemc_rad);
  tcm->setdphicut(0.15);
  tcm->setdzcut(10);
  tcm->setTrackPtLowCut(0.2);
  tcm->setEmcalELowCut(0.1);
  tcm->setnTpcClusters(20);
  tcm->setTrackQuality(1000);
  tcm->setRawClusContEMName("CLUSTERINFO_CEMC");
  tcm->setRawTowerGeomContName("TOWERGEOM_CEMCv3");
  se->registerSubsystem(tcm);

  TString photonconv_kfp_likesign_outfile = theOutfile + "_photonconv_kfp_likesign.root";
  std::string photonconv_kfp_likesign_string(photonconv_kfp_likesign_outfile.Data());
  KFPReco("PhotonConvKFPReco_likesign", "[gamma -> e^+ e^+]cc", photonconv_kfp_likesign_string, "MySvtxTrackMap", "PhotonConv_likesign");

  TString photonconv_kfp_unlikesign_outfile = theOutfile + "_photonconv_kfp_unlikesign.root";
  std::string photonconv_kfp_unlikesign_string(photonconv_kfp_unlikesign_outfile.Data());
  KFPReco("PhotonConvKFPReco_unlikesign", "gamma -> e^+ e^-", photonconv_kfp_unlikesign_string, "MySvtxTrackMap", "PhotonConv_unlikesign");

  TString track2calo_unlikesign_outfile = theOutfile + "_track2calo_unlikesign.root";
  std::string track2calo_unlikesign_string(track2calo_unlikesign_outfile.Data());
  TrackToCalo *ttc_unlikesign = new TrackToCalo("Tracks_And_Calo", track2calo_unlikesign_string);
  ttc_unlikesign->EMcalRadiusUser(doEMcalRadiusCorr);
  ttc_unlikesign->setEMcalRadius(new_cemc_rad);
  ttc_unlikesign->doLikesign(false);
  ttc_unlikesign->setKFPtrackMapName("PhotonConv_unlikesign_SvtxTrackMap");
  ttc_unlikesign->setKFPContName("PhotonConv_unlikesign_KFParticle_Container");
  ttc_unlikesign->anaTrkrInfo(false); // general track QA
  ttc_unlikesign->anaCaloInfo(false); // general calo QA
  ttc_unlikesign->doTrkrCaloMatching(false); // SvtxTrack match with calo
  ttc_unlikesign->doTrkrCaloMatching_KFP(true); // KFP selected trck match with calo
  ttc_unlikesign->setTrackPtLowCut(0.2);
  ttc_unlikesign->setEmcalELowCut(0.1);
  ttc_unlikesign->setnTpcClusters(20);
  ttc_unlikesign->setTrackQuality(1000);
  ttc_unlikesign->setRawClusContEMName("CLUSTERINFO_CEMC");
  ttc_unlikesign->setRawTowerGeomContName("TOWERGEOM_CEMCv3");
  se->registerSubsystem(ttc_unlikesign); 

  TString track2calo_likesign_outfile = theOutfile + "_track2calo_likesign.root";
  std::string track2calo_likesign_string(track2calo_likesign_outfile.Data());
  TrackToCalo *ttc_likesign = new TrackToCalo("Tracks_And_Calo", track2calo_likesign_string);
  ttc_likesign->EMcalRadiusUser(doEMcalRadiusCorr);
  ttc_likesign->setEMcalRadius(new_cemc_rad);
  ttc_likesign->doLikesign(true);
  ttc_likesign->setKFPtrackMapName("PhotonConv_likesign_SvtxTrackMap");
  ttc_likesign->setKFPContName("PhotonConv_likesign_KFParticle_Container");
  ttc_likesign->anaTrkrInfo(false); // general track QA
  ttc_likesign->anaCaloInfo(false); // general calo QA
  ttc_likesign->doTrkrCaloMatching(false); // SvtxTrack match with calo
  ttc_likesign->doTrkrCaloMatching_KFP(true); // KFP selected trck match with calo
  ttc_likesign->setTrackPtLowCut(0.2);
  ttc_likesign->setEmcalELowCut(0.1);
  ttc_likesign->setnTpcClusters(20);
  ttc_likesign->setTrackQuality(1000);
  ttc_likesign->setRawClusContEMName("CLUSTERINFO_CEMC");
  ttc_likesign->setRawTowerGeomContName("TOWERGEOM_CEMCv3");
  se->registerSubsystem(ttc_likesign); 

  se->skip(stepsize*index);
  se->run(nEvents);
  se->End();
  se->PrintTimer();
  std::cout << "CDB Files used:" << std::endl;
  CDBInterface::instance()->Print();

  ifstream file_photonconv_kfp_likesign(photonconv_kfp_likesign_string.c_str(), ios::binary | ios::ate);
  if (file_photonconv_kfp_likesign.good() && (file_photonconv_kfp_likesign.tellg() > 100))
  {
    string outputDirMove = outdir + "/Reconstructed/" + to_string(runnumber) + "/";
    string makeDirectoryMove = "mkdir -p " + outputDirMove;
    system(makeDirectoryMove.c_str());
    string moveOutput = "mv " + photonconv_kfp_likesign_string + " " + outputDirMove;
    std::cout << "moveOutput: " << moveOutput << std::endl;
    system(moveOutput.c_str());
  }

  ifstream file_photonconv_kfp_unlikesign(photonconv_kfp_unlikesign_string.c_str(), ios::binary | ios::ate);
  if (file_photonconv_kfp_unlikesign.good() && (file_photonconv_kfp_unlikesign.tellg() > 100))
  {
    string outputDirMove = outdir + "/Reconstructed/" + to_string(runnumber) + "/";
    string makeDirectoryMove = "mkdir -p " + outputDirMove;
    system(makeDirectoryMove.c_str());
    string moveOutput = "mv " + photonconv_kfp_unlikesign_string + " " + outputDirMove;
    std::cout << "moveOutput: " << moveOutput << std::endl;
    system(moveOutput.c_str());
  }

  ifstream file_track2calo_unlikesign(track2calo_unlikesign_string.c_str(), ios::binary | ios::ate);
  if (file_track2calo_unlikesign.good() && (file_track2calo_unlikesign.tellg() > 100))
  {
    string outputDirMove = outdir + "/Reconstructed/" + to_string(runnumber) + "/";
    string makeDirectoryMove = "mkdir -p " + outputDirMove;
    system(makeDirectoryMove.c_str());
    string moveOutput = "mv " + track2calo_unlikesign_string + " " + outputDirMove;
    std::cout << "moveOutput: " << moveOutput << std::endl;
    system(moveOutput.c_str());
  }

  ifstream file_track2calo_likesign(track2calo_likesign_string.c_str(), ios::binary | ios::ate);
  if (file_track2calo_likesign.good() && (file_track2calo_likesign.tellg() > 100))
  {
    string outputDirMove = outdir + "/Reconstructed/" + to_string(runnumber) + "/";
    string makeDirectoryMove = "mkdir -p " + outputDirMove;
    system(makeDirectoryMove.c_str());
    string moveOutput = "mv " + track2calo_likesign_string + " " + outputDirMove;
    std::cout << "moveOutput: " << moveOutput << std::endl;
    system(moveOutput.c_str());
  }

  delete se;
  std::cout << "Finished" << std::endl;
  gSystem->Exit(0);
}

void KFPReco(std::string module_name = "KFPReco", std::string decaydescriptor = "K_S0 -> pi^+ pi^-", std::string outfile = "KFP.root", std::string trackmapName = "SvtxTrackMap", std::string containerName = "KFParticle")
{
  auto se = Fun4AllServer::instance();
  //KFParticle setup
  KFParticle_sPHENIX *kfparticle = new KFParticle_sPHENIX(module_name);
  kfparticle->Verbosity(0);
  kfparticle->setDecayDescriptor(decaydescriptor);

  kfparticle->setTrackMapNodeName(trackmapName);
  kfparticle->setContainerName(containerName);

  //Basic node selection and configuration
  kfparticle->magFieldFile("FIELDMAP_TRACKING");
  kfparticle->getAllPVInfo(false);
  kfparticle->allowZeroMassTracks(true);
  kfparticle->getDetectorInfo(true);
  kfparticle->useFakePrimaryVertex(false);
  kfparticle->saveDST();

  kfparticle->constrainToPrimaryVertex(true);
  kfparticle->setMotherIPchi2(FLT_MAX);
  kfparticle->setFlightDistancechi2(-1.);
  kfparticle->setMinDIRA(-1.1);
  kfparticle->setDecayLengthRange(0., FLT_MAX);
  kfparticle->setDecayTimeRange(-1*FLT_MAX, FLT_MAX);

  //Track parameters
  kfparticle->setMinMVTXhits(0);
  //kfparticle->setMinTPChits(20);
  kfparticle->setMinTPChits(0);
  kfparticle->setMinimumTrackPT(0.2);
  kfparticle->setMaximumTrackPTchi2(FLT_MAX);
  kfparticle->setMinimumTrackIPchi2(-1.);
  kfparticle->setMinimumTrackIP(-1.);
  //kfparticle->setMaximumTrackchi2nDOF(100.);
  kfparticle->setMaximumTrackchi2nDOF(FLT_MAX);

  //Vertex parameters
  //kfparticle->setMaximumVertexchi2nDOF(50);
  kfparticle->setMaximumVertexchi2nDOF(FLT_MAX);
  //kfparticle->setMaximumDaughterDCA(1.);
  kfparticle->setMaximumDaughterDCA(FLT_MAX);

  //Parent parameters
  kfparticle->setMotherPT(0);
  kfparticle->setMinimumMass(-1);
  kfparticle->setMaximumMass(10);
  //kfparticle->setMaximumMotherVertexVolume(0.1);
  kfparticle->setMaximumMotherVertexVolume(FLT_MAX);

  kfparticle->setOutputName(outfile);

  se->registerSubsystem(kfparticle);
}
