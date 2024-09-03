#include <fun4all/Fun4AllServer.h>

#pragma GCC diagnostic push

#pragma GCC diagnostic ignored "-Wundefined-internal"

#include <kfparticle_sphenix/KFParticle_sPHENIX.h>

#pragma GCC diagnostic pop

#include <phhepmc/Fun4AllHepMCInputManager.h>
#include <phhepmc/Fun4AllHepMCPileupInputManager.h>
#include <phhepmc/HepMCFlowAfterBurner.h>
#include <phhepmc/PHHepMCGenHelper.h>

#include <GlobalVariables.C>

#include <G4_Mbd.C>
#include <G4_BlackHole.C>
#include <G4_CEmc_Albedo.C>
#include <G4_CEmc_Spacal.C>
#include <G4_EPD.C>
#include <G4_HcalIn_ref.C>
#include <G4_HcalOut_ref.C>
#include <G4_BeamLine.C>
#include <G4_Magnet.C>
#include <G4_PSTOF.C>
#include <G4_Pipe.C>
#include <G4_PlugDoor.C>
#include <G4_TrkrSimulation.C>
#include <G4_User.C>
#include <G4_World.C>
#include <G4_ZDC.C>

#include <g4main/PHG4Reco.h>
#include <g4main/PHG4TruthSubsystem.h>

R__LOAD_LIBRARY(libkfparticle_sphenix.so)

namespace HeavyFlavorReco
{
  // https://wiki.bnl.gov/sPHENIX/index.php/KFParticle
  string decayDescriptor = "gamma -> e^+ e^-";  //See twiki on how to set this
  //string decayDescriptor = "[D0 -> K^- pi^+]cc";  //See twiki on how to set this
  string reconstructionName = "PhotonConvKFPReco";         //Used for naming output folder, file and node
  string trackmapName = "MySvtxTrackMap";
  //Used for naming container in DST
  string containerName = "PhotonConv";
  string outputRecoFile;
  string outputHFEffFile;
  string outputTrackingEvalFile;
  bool useMyTrackMap = true;  //Alternative Track Map
  bool useContainer = true; // Save Container: svtxtrack and KFParticle
  bool runTruthTrigger = false;  //Decay Finder
  bool runTrackEff = false;  //HF track efficiency
  bool saveEvalFile = false;  //Official Eval root file
  bool getTruthInfo = false;      //Add truth matching to output file
  bool getCaloInfo = false;
  bool runTracking = true;  //Run tracking on DSTs
  bool buildTruthTable = true;
  bool runQA = false;        //Run QA, needs set up
  int VERBOSITY = 0;

  //! apply reference sPHENIX nominal beam parameter with 2mrad crossing as defined in sPH-TRG-2022-001 and past RHIC experience
  //! \param[in] HepMCGen any HepMC generator, e.g. Fun4AllHepMCInputManager, Fun4AllHepMCPileupInputManager, PHPythia8, PHPythia6, ReadEICFiles
  //! \param[in] collision_type select the beam configuration with Input::BeamConfiguration
  void ApplysPHENIXBeamParameter(PHHepMCGenHelper *HepMCGen, const Input::BeamConfiguration & beam_config)
  {
    if (HepMCGen == nullptr)
    {
      std::cout << "ApplysPHENIXBeamParameter(): Fatal Error - null input pointer HepMCGen" << std::endl;
      exit(1);
    }
    HepMCGen->set_beam_direction_theta_phi(1e-3, 0, M_PI - 1e-3, 0);  //2mrad x-ing of sPHENIX per sPH-TRG-2022-001

    switch (beam_config)
    {
    case Input::AA_COLLISION:
      // heavy ion mode

      HepMCGen->set_vertex_distribution_width(
          100e-4,         // approximation from past STAR/Run16 AuAu data
          100e-4,         // approximation from past STAR/Run16 AuAu data
          7,              // sPH-TRG-2022-001. Fig B.2
          20 / 29.9792);  // 20cm collision length / speed of light in cm/ns

      break;
    case Input::pA_COLLISION:

      // pA mode

      HepMCGen->set_vertex_distribution_width(
          100e-4,         // set to be similar to AA
          100e-4,         // set to be similar to AA
          8,              // sPH-TRG-2022-001. Fig B.4
          20 / 29.9792);  // 20cm collision length / speed of light in cm/ns

      break;
    case Input::pp_COLLISION:

      // pp mode

      HepMCGen->set_vertex_distribution_width(
          120e-4,         // approximation from past PHENIX data
          120e-4,         // approximation from past PHENIX data
          10,              // sPH-TRG-2022-001. Fig B.3
          20 / 29.9792);  // 20cm collision length / speed of light in cm/ns

      break;
    default:
      std::cout <<"ApplysPHENIXBeamParameter: invalid beam_config = "<<beam_config<<std::endl;

      exit(1);

    }

    HepMCGen->set_vertex_distribution_function(
        PHHepMCGenHelper::Gaus,
        PHHepMCGenHelper::Gaus,
        PHHepMCGenHelper::Gaus,
        PHHepMCGenHelper::Gaus);
  }

  //! apply sPHENIX nominal beam parameter according to the beam collision setting of Input::IS_PP_COLLISION
  //! \param[in] HepMCGen any HepMC generator, e.g. Fun4AllHepMCInputManager, Fun4AllHepMCPileupInputManager, PHPythia8, PHPythia6, ReadEICFiles
  void ApplysPHENIXBeamParameter(PHHepMCGenHelper *HepMCGen)
  {
    ApplysPHENIXBeamParameter(HepMCGen, Input::BEAM_CONFIGURATION);
  }

  void G4Init()
  {
    PipeInit();
    MvtxInit();
    InttInit();
    TPCInit();
    MicromegasInit();

    G4MAGNET::magfield_rescale = 1.;
    MagnetInit();
    MagnetFieldInit();
  }

  int G4Setup()
  {
  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer *se = Fun4AllServer::instance();

  PHG4Reco *g4Reco = new PHG4Reco();
  g4Reco->set_rapidity_coverage(1.1);  // according to drawings
  WorldInit(g4Reco);
  //PYTHIA 6
  if (G4P6DECAYER::decayType != EDecayType::kAll)
  {
    g4Reco->set_force_decay(G4P6DECAYER::decayType);
  }
  //EvtGen 
  g4Reco->CustomizeEvtGenDecay(EVTGENDECAYER::DecayFile); 

  double fieldstrength;
  istringstream stringline(G4MAGNET::magfield);
  stringline >> fieldstrength;
  if (stringline.fail())
  {  // conversion to double fails -> we have a string

    if (G4MAGNET::magfield.find("sphenix3dbigmapxyz") != string::npos ||
        G4MAGNET::magfield == "CDB")
    {
      g4Reco->set_field_map(G4MAGNET::magfield, PHFieldConfig::Field3DCartesian);
    }
    else
    {
      g4Reco->set_field_map(G4MAGNET::magfield, PHFieldConfig::kField2D);
    }
  }
  else
  {
    g4Reco->set_field(fieldstrength);  // use const soleniodal field
  }
  g4Reco->set_field_rescale(G4MAGNET::magfield_rescale);

// the radius is an older protection against overlaps, it is not
// clear how well this works nowadays but it doesn't hurt either
  double radius = 0.;

  radius = Pipe(g4Reco, radius);
  radius = Mvtx(g4Reco, radius);
  radius = Intt(g4Reco, radius);
  radius = TPC(g4Reco, radius);
  Micromegas(g4Reco);

  PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);

  // finally adjust the world size in case the default is too small
  WorldSize(g4Reco, radius);

  se->registerSubsystem(g4Reco);
  return 0;
  }

};  // namespace HeavyFlavorReco

using namespace HeavyFlavorReco;

void myHeavyFlavorReco()
{
  Fun4AllServer *se = Fun4AllServer::instance();

  KFParticle_sPHENIX *kfparticle = new KFParticle_sPHENIX(reconstructionName);
  kfparticle->Verbosity(1);

  kfparticle->setDecayDescriptor(decayDescriptor);

  if (runTrackEff) kfparticle->setTrackMapNodeName("HFSelected_SvtxTrackMap");

  kfparticle->doTruthMatching(getTruthInfo);
  kfparticle->getDetectorInfo(false);
  kfparticle->getCaloInfo(getCaloInfo);
  kfparticle->getAllPVInfo(false);
  kfparticle->allowZeroMassTracks(true);
  kfparticle->saveDST(false);
  kfparticle->saveParticleContainer(false);

  bool fixToPV = false;
  bool useFakePV = false;

  if (useFakePV)
  {
    fixToPV = false;  //Constraining to a fake PV results in some gibberish variables
    kfparticle->useFakePrimaryVertex(true);
  }

  if (fixToPV)
  {
    kfparticle->constrainToPrimaryVertex(true);
    kfparticle->setMotherIPchi2(3);
    kfparticle->setFlightDistancechi2(-1.);
    kfparticle->setMinDIRA(0.90);
    kfparticle->setDecayLengthRange(-1*FLT_MAX, FLT_MAX);
  }

  //Track parameters
  kfparticle->setMinimumTrackPT(0.2);
  kfparticle->setMinimumTrackIPchi2(-1.);
  kfparticle->setMinimumTrackIP(-1.);
  kfparticle->setMaximumTrackchi2nDOF(100);

  //Vertex parameters
  kfparticle->setMaximumVertexchi2nDOF(100);
  kfparticle->setMaximumDaughterDCA(0.1);

  //Parent parameters
  kfparticle->setMotherPT(0);
  kfparticle->setMinimumMass(0);
  kfparticle->setMaximumMass(0.5);
  kfparticle->setMaximumMotherVertexVolume(0.03);

  kfparticle->setOutputName(outputRecoFile);

  se->registerSubsystem(kfparticle);
}

void PhotonConvKFPReco()
{
  Fun4AllServer *se = Fun4AllServer::instance();

  KFParticle_sPHENIX *kfparticle = new KFParticle_sPHENIX(reconstructionName);
  kfparticle->Verbosity(1);

  if (useMyTrackMap) kfparticle->setTrackMapNodeName(trackmapName);
  if (useContainer) kfparticle->setContainerName(containerName);

  kfparticle->setDecayDescriptor(decayDescriptor);

  if (runTrackEff) kfparticle->setTrackMapNodeName("HFSelected_SvtxTrackMap");

  kfparticle->doTruthMatching(getTruthInfo);
  kfparticle->getDetectorInfo(false);
  kfparticle->getCaloInfo(getCaloInfo);
  kfparticle->getAllPVInfo(false);
  kfparticle->allowZeroMassTracks(true);
  kfparticle->saveDST(true);
  //kfparticle->saveParticleContainer(false);

  bool fixToPV = false;
  bool useFakePV = true;

  if (useFakePV)
  {
    fixToPV = false;  //Constraining to a fake PV results in some gibberish variables
    kfparticle->useFakePrimaryVertex(true);
  }

  if (fixToPV)
  {
    kfparticle->constrainToPrimaryVertex(true);
    kfparticle->setMotherIPchi2(FLT_MAX);
    kfparticle->setFlightDistancechi2(-1.);
    kfparticle->setMinDIRA(-1.1);
    kfparticle->setDecayLengthRange(-1*FLT_MAX, FLT_MAX);
  }

  //Track parameters
  kfparticle->setMinimumTrackPT(0.);
  kfparticle->setMinimumTrackIPchi2(-1.);
  kfparticle->setMinimumTrackIP(-1.);
  kfparticle->setMaximumTrackchi2nDOF(FLT_MAX);

  //Vertex parameters
  kfparticle->setMaximumVertexchi2nDOF(FLT_MAX);
  kfparticle->setMaximumDaughterDCA(10);

  //Parent parameters
  kfparticle->setMotherPT(0);
  kfparticle->setMinimumMass(-1.0);
  kfparticle->setMaximumMass(10.0);
  kfparticle->setMaximumMotherVertexVolume(0.9);
  kfparticle->setMotherIPchi2(FLT_MAX);
  kfparticle->setFlightDistancechi2(-1.);
  kfparticle->setMinDIRA(-1.1);
  kfparticle->setDecayLengthRange(-1*FLT_MAX, FLT_MAX);

  kfparticle->setOutputName(outputRecoFile);

  se->registerSubsystem(kfparticle);
}

void myDemoReco()
{
  Fun4AllServer *se = Fun4AllServer::instance();

  KFParticle_sPHENIX *kfparticle = new KFParticle_sPHENIX(reconstructionName);
  kfparticle->Verbosity(1);

  kfparticle->setDecayDescriptor("B_s0 -> {J/psi -> e^+ e^-} {K_S0 -> pi^+ pi^-}");
  kfparticle->setTrackMapNodeName("HFSelected_SvtxTrackMap");

  kfparticle->doTruthMatching(true);
  kfparticle->getCaloInfo(false);
  kfparticle->allowZeroMassTracks(true);
  kfparticle->constrainToPrimaryVertex(true);
  kfparticle->setOutputName("kfparticle_demo.root");

  //Track parameters
  kfparticle->setMinimumTrackPT(0.);
  kfparticle->setMinimumTrackIPchi2(-1.);
  kfparticle->setMinimumTrackIP(-1.);
  kfparticle->setMaximumTrackchi2nDOF(90.);

  //Vertex parameters
  kfparticle->setMaximumVertexchi2nDOF(90);
  kfparticle->setMaximumDaughterDCA(0.7);

  //Parent parameters
  kfparticle->setMotherPT(0);
  kfparticle->setMinimumMass(2.0);
  kfparticle->setMaximumMass(6.0);
  kfparticle->setMaximumMotherVertexVolume(0.9);
  kfparticle->setMotherIPchi2(FLT_MAX);
  kfparticle->setFlightDistancechi2(-1.);
  kfparticle->setMinDIRA(-1.1);
  kfparticle->setDecayLengthRange(-1*FLT_MAX, FLT_MAX);

  //Intermediate parameters
  std::vector<std::pair<float, float>> intermediate_mass_range;
  intermediate_mass_range.push_back(make_pair(0.8, 3.5));
  intermediate_mass_range.push_back(make_pair(0.4, 0.6));
  kfparticle->setIntermediateMassRange(intermediate_mass_range);

  std::vector<float> intermediate_min_pt = {0.0, 0.0};
  kfparticle->setIntermediateMinPT(intermediate_min_pt);

  std::vector<std::pair<float, float>> intermediate_IP_range;
  intermediate_IP_range.push_back(make_pair(-1., 5.));
  intermediate_IP_range.push_back(make_pair(-1., 5.));
  kfparticle->setIntermediateIPRange(intermediate_IP_range);

  std::vector<std::pair<float, float>> intermediate_IPchi2_range;
  intermediate_IPchi2_range.push_back(make_pair(0., 400.));
  intermediate_IPchi2_range.push_back(make_pair(0., 400.));
  kfparticle->setIntermediateIPchi2Range(intermediate_IPchi2_range);

  std::vector<float> intermediate_min_dira = {-1.1, -1.1};
  kfparticle->setIntermediateMinDIRA(intermediate_min_dira);

  std::vector<float> intermediate_min_FDchi2 = {-1., -1.};
  kfparticle->setIntermediateMinFDchi2(intermediate_min_FDchi2);

  std::vector<float> intermediate_max_vertex_vol = {1.1, 0.9};
  kfparticle->setIntermediateMaxVertexVolume(intermediate_max_vertex_vol);

  se->registerSubsystem(kfparticle);
}
