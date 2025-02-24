/*!
 *  \file   TrackToCalo.cc
 *  \brief  Track To Calo, output root file
 *  \author Xudong Yu <xyu3@bnl.gov>
 */
#include "TrackToCalo.h"

#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoDefs.h>

#include <ffarawobjects/Gl1Packet.h>
#include <ffaobjects/EventHeaderv1.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/SvtxVertexMap.h>
#include <globalvertex/SvtxVertex.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterCrossingAssocv1.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackAnalysisUtils.h>
#include <trackreco/ActsPropagator.h>

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Surfaces/CylinderSurface.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>

#include <KFVertex.h>

#include <CLHEP/Vector/ThreeVector.h>
#include <math.h>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include <kfparticle_sphenix/KFParticle_Tools.h>
KFParticle_Tools kf_tools;

namespace
{
  //! get cluster keys from a given track
  std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack* track)
  {
    std::vector<TrkrDefs::cluskey> out;
    for (const auto& seed : {track->get_silicon_seed(), track->get_tpc_seed()})
    {
      if (seed)
      {
        std::copy(seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter(out));
      }
    }
    return out;
  }

  /// return number of clusters of a given type that belong to a tracks
  template <int type>
  int count_clusters(const std::vector<TrkrDefs::cluskey>& keys)
  {
    return std::count_if(keys.begin(), keys.end(),
                         [](const TrkrDefs::cluskey& key)
                         { return TrkrDefs::getTrkrId(key) == type; });
  }
}

//____________________________________________________________________________..
TrackToCalo::TrackToCalo(const std::string &name, const std::string &file):
 SubsysReco(name),
 _outfilename(file),
 _outfile(nullptr),
 _tree(nullptr),
 _tree_KFP(nullptr)
{
  std::cout << "TrackToCalo::TrackToCalo(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
TrackToCalo::~TrackToCalo()
{
  std::cout << "TrackToCalo::~TrackToCalo() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int TrackToCalo::Init(PHCompositeNode *topNode)
{
  std::cout << topNode << std::endl;
  std::cout << "TrackToCalo::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  delete _outfile;
  _outfile = new TFile(_outfilename.c_str(), "RECREATE");

  if (m_doTrkrCaloMatching)
  {
    createBranches();
  }

  if (m_doTrkrCaloMatching_KFP)
  {
    createBranches_KFP();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void TrackToCalo::createBranches()
{
  delete _tree;
  _tree = new TTree("tree", "A tree with track/calo info");
  _tree->Branch("_runNumber", &_runNumber);
  _tree->Branch("_eventNumber", &_eventNumber);
  _tree->Branch("_vertex_id", &_vertex_id);
  _tree->Branch("_vertex_crossing", &_vertex_crossing);
  _tree->Branch("_vertex_ntracks", &_vertex_ntracks);
  _tree->Branch("_vertex_x", &_vertex_x);
  _tree->Branch("_vertex_y", &_vertex_y);
  _tree->Branch("_vertex_z", &_vertex_z);
  _tree->Branch("_cluster_x", &_cluster_x);
  _tree->Branch("_cluster_y", &_cluster_y);
  _tree->Branch("_cluster_z", &_cluster_z);
  _tree->Branch("_track_id", &_track_id);
  _tree->Branch("_track_bc", &_track_bc);
  _tree->Branch("_track_phi", &_track_phi);
  _tree->Branch("_track_eta", &_track_eta);
  _tree->Branch("_track_pcax", &_track_pcax);
  _tree->Branch("_track_pcay", &_track_pcay);
  _tree->Branch("_track_pcaz", &_track_pcaz);
  _tree->Branch("_track_crossing", &_track_crossing);
  _tree->Branch("_track_vx", &_track_vx);
  _tree->Branch("_track_vy", &_track_vy);
  _tree->Branch("_track_vz", &_track_vz);
  _tree->Branch("_track_quality", &_track_quality);
  _tree->Branch("_track_dcaxy", &_track_dcaxy);
  _tree->Branch("_track_dcaz", &_track_dcaz);
  _tree->Branch("_track_nc_mvtx", &_track_nc_mvtx);
  _tree->Branch("_track_nc_intt", &_track_nc_intt);
  _tree->Branch("_track_nc_tpc", &_track_nc_tpc);
  _tree->Branch("_track_nc_tpot", &_track_nc_tpot);
  _tree->Branch("_track_ptq", &_track_ptq);
  _tree->Branch("_track_px", &_track_px);
  _tree->Branch("_track_py", &_track_py);
  _tree->Branch("_track_pz", &_track_pz);
  _tree->Branch("_track_phi_origin", &_track_phi_origin);
  _tree->Branch("_track_eta_origin", &_track_eta_origin);
  _tree->Branch("_track_px_origin", &_track_px_origin);
  _tree->Branch("_track_py_origin", &_track_py_origin);
  _tree->Branch("_track_pz_origin", &_track_pz_origin);
  _tree->Branch("_track_x_origin", &_track_x_origin);
  _tree->Branch("_track_y_origin", &_track_y_origin);
  _tree->Branch("_track_z_origin", &_track_z_origin);
  _tree->Branch("_track_phi_emc", &_track_phi_emc);
  _tree->Branch("_track_eta_emc", &_track_eta_emc);
  _tree->Branch("_track_px_emc", &_track_px_emc);
  _tree->Branch("_track_py_emc", &_track_py_emc);
  _tree->Branch("_track_pz_emc", &_track_pz_emc);
  _tree->Branch("_track_x_emc", &_track_x_emc);
  _tree->Branch("_track_y_emc", &_track_y_emc);
  _tree->Branch("_track_z_emc", &_track_z_emc);
  _tree->Branch("_track_phi_ihc", &_track_phi_ihc);
  _tree->Branch("_track_eta_ihc", &_track_eta_ihc);
  _tree->Branch("_track_px_ihc", &_track_px_ihc);
  _tree->Branch("_track_py_ihc", &_track_py_ihc);
  _tree->Branch("_track_pz_ihc", &_track_pz_ihc);
  _tree->Branch("_track_x_ihc", &_track_x_ihc);
  _tree->Branch("_track_y_ihc", &_track_y_ihc);
  _tree->Branch("_track_z_ihc", &_track_z_ihc);
  _tree->Branch("_track_phi_ohc", &_track_phi_ohc);
  _tree->Branch("_track_eta_ohc", &_track_eta_ohc);
  _tree->Branch("_track_px_ohc", &_track_px_ohc);
  _tree->Branch("_track_py_ohc", &_track_py_ohc);
  _tree->Branch("_track_pz_ohc", &_track_pz_ohc);
  _tree->Branch("_track_x_ohc", &_track_x_ohc);
  _tree->Branch("_track_y_ohc", &_track_y_ohc);
  _tree->Branch("_track_z_ohc", &_track_z_ohc);
  _tree->Branch("_trClus_track_id", &_trClus_track_id);
  _tree->Branch("_trClus_type", &_trClus_type);
  _tree->Branch("_trClus_x", &_trClus_x);
  _tree->Branch("_trClus_y", &_trClus_y);
  _tree->Branch("_trClus_z", &_trClus_z);
  _tree->Branch("_emcal_id", &_emcal_id);
  _tree->Branch("_emcal_phi", &_emcal_phi);
  _tree->Branch("_emcal_eta", &_emcal_eta);
  _tree->Branch("_emcal_x", &_emcal_x);
  _tree->Branch("_emcal_y", &_emcal_y);
  _tree->Branch("_emcal_z", &_emcal_z);
  _tree->Branch("_emcal_e", &_emcal_e);
  _tree->Branch("_emcal_ecore", &_emcal_ecore);
  _tree->Branch("_emcal_chi2", &_emcal_chi2);
  _tree->Branch("_emcal_prob", &_emcal_prob);
  _tree->Branch("_emcal_tower_cluster_id", &_emcal_tower_cluster_id);
  _tree->Branch("_emcal_tower_e", &_emcal_tower_e);
  _tree->Branch("_emcal_tower_phi", &_emcal_tower_phi);
  _tree->Branch("_emcal_tower_eta", &_emcal_tower_eta);
  _tree->Branch("_emcal_tower_status", &_emcal_tower_status);
  _tree->Branch("_hcal_id", &_hcal_id);
  _tree->Branch("_hcal_phi", &_hcal_phi);
  _tree->Branch("_hcal_eta", &_hcal_eta);
  _tree->Branch("_hcal_x", &_hcal_x);
  _tree->Branch("_hcal_y", &_hcal_y);
  _tree->Branch("_hcal_z", &_hcal_z);
  _tree->Branch("_hcal_e", &_hcal_e);
  _tree->Branch("_hcal_tower_cluster_id", &_hcal_tower_cluster_id);
  _tree->Branch("_hcal_tower_e", &_hcal_tower_e);
  _tree->Branch("_hcal_tower_phi", &_hcal_tower_phi);
  _tree->Branch("_hcal_tower_eta", &_hcal_tower_eta);
  _tree->Branch("_hcal_tower_status", &_hcal_tower_status);
  _tree->Branch("_hcal_tower_io", &_hcal_tower_io);
  _tree->Branch("_mbd_x", &_mbd_x);
  _tree->Branch("_mbd_y", &_mbd_y);
  _tree->Branch("_mbd_z", &_mbd_z);
  _tree->Branch("_triggers", &_triggers);
  _tree->Branch("_ntracks", &_ntracks);
}

void TrackToCalo::createBranches_KFP()
{
  delete _tree_KFP;
  _tree_KFP = new TTree("tree_KFP", "A tree with track/calo info after KFParticle");
  _tree_KFP->Branch("_runNumber", &_runNumber);
  _tree_KFP->Branch("_eventNumber", &_eventNumber);
  _tree_KFP->Branch("_numCan", &_numCan);
  _tree_KFP->Branch("_gamma_mass", &_gamma_mass);
  _tree_KFP->Branch("_gamma_massErr", &_gamma_massErr);
  _tree_KFP->Branch("_gamma_charge", &_gamma_charge);
  _tree_KFP->Branch("_gamma_x", &_gamma_x);
  _tree_KFP->Branch("_gamma_y", &_gamma_y);
  _tree_KFP->Branch("_gamma_z", &_gamma_z);
  _tree_KFP->Branch("_gamma_px", &_gamma_px);
  _tree_KFP->Branch("_gamma_py", &_gamma_py);
  _tree_KFP->Branch("_gamma_pz", &_gamma_pz);
  _tree_KFP->Branch("_gamma_pE", &_gamma_pE);
  _tree_KFP->Branch("_gamma_pT", &_gamma_pT);
  _tree_KFP->Branch("_gamma_pTErr", &_gamma_pTErr);
  _tree_KFP->Branch("_gamma_p", &_gamma_p);
  _tree_KFP->Branch("_gamma_pErr", &_gamma_pErr);
  _tree_KFP->Branch("_gamma_pseudorapidity", &_gamma_pseudorapidity);
  _tree_KFP->Branch("_gamma_rapidity", &_gamma_rapidity);
  _tree_KFP->Branch("_gamma_theta", &_gamma_theta);
  _tree_KFP->Branch("_gamma_phi", &_gamma_phi);
  _tree_KFP->Branch("_gamma_chi2", &_gamma_chi2);
  _tree_KFP->Branch("_gamma_nDoF", &_gamma_nDoF);
  _tree_KFP->Branch("_gamma_vertex_volume", &_gamma_vertex_volume);
  _tree_KFP->Branch("_gamma_SV_chi2_per_nDoF", &_gamma_SV_chi2_per_nDoF);

  _tree_KFP->Branch("_ep_mass", &_ep_mass);
  _tree_KFP->Branch("_ep_x", &_ep_x);
  _tree_KFP->Branch("_ep_x_raw", &_ep_x_raw);
  _tree_KFP->Branch("_ep_y", &_ep_y);
  _tree_KFP->Branch("_ep_y_raw", &_ep_y_raw);
  _tree_KFP->Branch("_ep_z", &_ep_z);
  _tree_KFP->Branch("_ep_z_raw", &_ep_z_raw);
  _tree_KFP->Branch("_ep_px", &_ep_px);
  _tree_KFP->Branch("_ep_px_raw", &_ep_px_raw);
  _tree_KFP->Branch("_ep_py", &_ep_py);
  _tree_KFP->Branch("_ep_py_raw", &_ep_py_raw);
  _tree_KFP->Branch("_ep_pz", &_ep_pz);
  _tree_KFP->Branch("_ep_pz_raw", &_ep_pz_raw);
  _tree_KFP->Branch("_ep_pE", &_ep_pE);
  _tree_KFP->Branch("_ep_pE_unmoved", &_ep_pE_unmoved);
  _tree_KFP->Branch("_ep_pT", &_ep_pT);
  _tree_KFP->Branch("_ep_pTErr", &_ep_pTErr);
  _tree_KFP->Branch("_ep_pT_raw", &_ep_pT_raw);
  _tree_KFP->Branch("_ep_pT_unmoved", &_ep_pT_unmoved);
  _tree_KFP->Branch("_ep_p", &_ep_p);
  _tree_KFP->Branch("_ep_pErr", &_ep_pErr);
  _tree_KFP->Branch("_ep_p_raw", &_ep_p_raw);
  _tree_KFP->Branch("_ep_p_unmoved", &_ep_p_unmoved);
  _tree_KFP->Branch("_ep_pseudorapidity", &_ep_pseudorapidity);
  _tree_KFP->Branch("_ep_pseudorapidity_raw", &_ep_pseudorapidity_raw);
  _tree_KFP->Branch("_ep_rapidity", &_ep_rapidity);
  _tree_KFP->Branch("_ep_theta", &_ep_theta);
  _tree_KFP->Branch("_ep_phi", &_ep_phi);
  _tree_KFP->Branch("_ep_phi_raw", &_ep_phi_raw);
  _tree_KFP->Branch("_ep_chi2", &_ep_chi2);
  _tree_KFP->Branch("_ep_chi2_raw", &_ep_chi2_raw);
  _tree_KFP->Branch("_ep_nDoF", &_ep_nDoF);
  _tree_KFP->Branch("_ep_nDoF_raw", &_ep_nDoF_raw);
  _tree_KFP->Branch("_ep_crossing", &_ep_crossing);
  _tree_KFP->Branch("_ep_nmvtx", &_ep_nmvtx);
  _tree_KFP->Branch("_ep_nintt", &_ep_nintt);
  _tree_KFP->Branch("_ep_ntpc", &_ep_ntpc);
  _tree_KFP->Branch("_ep_ntpot", &_ep_ntpot);
  _tree_KFP->Branch("_ep_clus_ican", &_ep_clus_ican);
  //_tree_KFP->Branch("_ep_clus_type", &_ep_clus_type);
  _tree_KFP->Branch("_ep_clus_x", &_ep_clus_x);
  _tree_KFP->Branch("_ep_clus_y", &_ep_clus_y);
  _tree_KFP->Branch("_ep_clus_z", &_ep_clus_z);
  _tree_KFP->Branch("_ep_phi_emc", &_ep_phi_emc);
  _tree_KFP->Branch("_ep_eta_emc", &_ep_eta_emc);
  _tree_KFP->Branch("_ep_px_emc", &_ep_px_emc);
  _tree_KFP->Branch("_ep_py_emc", &_ep_py_emc);
  _tree_KFP->Branch("_ep_pz_emc", &_ep_pz_emc);
  _tree_KFP->Branch("_ep_x_emc", &_ep_x_emc);
  _tree_KFP->Branch("_ep_y_emc", &_ep_y_emc);
  _tree_KFP->Branch("_ep_z_emc", &_ep_z_emc);
  if(m_doSimulation)
  {
    _tree_KFP->Branch("_ep_has_truthmatching", &_ep_has_truthmatching);
    _tree_KFP->Branch("_ep_true_id", &_ep_true_id);
    _tree_KFP->Branch("_ep_true_px", &_ep_true_px);
    _tree_KFP->Branch("_ep_true_py", &_ep_true_py);
    _tree_KFP->Branch("_ep_true_pz", &_ep_true_pz);
    _tree_KFP->Branch("_ep_true_vertex_x", &_ep_true_vertex_x);
    _tree_KFP->Branch("_ep_true_vertex_y", &_ep_true_vertex_y);
    _tree_KFP->Branch("_ep_true_vertex_z", &_ep_true_vertex_z);
    _tree_KFP->Branch("_ep_true_vertex_x_method2", &_ep_true_vertex_x_method2);
    _tree_KFP->Branch("_ep_true_vertex_y_method2", &_ep_true_vertex_y_method2);
    _tree_KFP->Branch("_ep_true_vertex_z_method2", &_ep_true_vertex_z_method2);
  }

  _tree_KFP->Branch("_em_mass", &_em_mass);
  _tree_KFP->Branch("_em_x", &_em_x);
  _tree_KFP->Branch("_em_x_raw", &_em_x_raw);
  _tree_KFP->Branch("_em_y", &_em_y);
  _tree_KFP->Branch("_em_y_raw", &_em_y_raw);
  _tree_KFP->Branch("_em_z", &_em_z);
  _tree_KFP->Branch("_em_z_raw", &_em_z_raw);
  _tree_KFP->Branch("_em_px", &_em_px);
  _tree_KFP->Branch("_em_px_raw", &_em_px_raw);
  _tree_KFP->Branch("_em_py", &_em_py);
  _tree_KFP->Branch("_em_py_raw", &_em_py_raw);
  _tree_KFP->Branch("_em_pz", &_em_pz);
  _tree_KFP->Branch("_em_pz_raw", &_em_pz_raw);
  _tree_KFP->Branch("_em_pE", &_em_pE);
  _tree_KFP->Branch("_em_pE_unmoved", &_em_pE_unmoved);
  _tree_KFP->Branch("_em_pT", &_em_pT);
  _tree_KFP->Branch("_em_pTErr", &_em_pTErr);
  _tree_KFP->Branch("_em_pT_raw", &_em_pT_raw);
  _tree_KFP->Branch("_em_pT_unmoved", &_em_pT_unmoved);
  _tree_KFP->Branch("_em_p", &_em_p);
  _tree_KFP->Branch("_em_pErr", &_em_pErr);
  _tree_KFP->Branch("_em_p_raw", &_em_p_raw);
  _tree_KFP->Branch("_em_p_unmoved", &_em_p_unmoved);
  _tree_KFP->Branch("_em_pseudorapidity", &_em_pseudorapidity);
  _tree_KFP->Branch("_em_pseudorapidity_raw", &_em_pseudorapidity_raw);
  _tree_KFP->Branch("_em_rapidity", &_em_rapidity);
  _tree_KFP->Branch("_em_theta", &_em_theta);
  _tree_KFP->Branch("_em_phi", &_em_phi);
  _tree_KFP->Branch("_em_phi_raw", &_em_phi_raw);
  _tree_KFP->Branch("_em_chi2", &_em_chi2);
  _tree_KFP->Branch("_em_chi2_raw", &_em_chi2_raw);
  _tree_KFP->Branch("_em_nDoF", &_em_nDoF);
  _tree_KFP->Branch("_em_nDoF_raw", &_em_nDoF_raw);
  _tree_KFP->Branch("_em_crossing", &_em_crossing);
  _tree_KFP->Branch("_em_nmvtx", &_em_nmvtx);
  _tree_KFP->Branch("_em_nintt", &_em_nintt);
  _tree_KFP->Branch("_em_ntpc", &_em_ntpc);
  _tree_KFP->Branch("_em_ntpot", &_em_ntpot);
  _tree_KFP->Branch("_em_clus_ican", &_em_clus_ican);
  //_tree_KFP->Branch("_em_clus_type", &_em_clus_type);
  _tree_KFP->Branch("_em_clus_x", &_em_clus_x);
  _tree_KFP->Branch("_em_clus_y", &_em_clus_y);
  _tree_KFP->Branch("_em_clus_z", &_em_clus_z);
  _tree_KFP->Branch("_em_phi_emc", &_em_phi_emc);
  _tree_KFP->Branch("_em_eta_emc", &_em_eta_emc);
  _tree_KFP->Branch("_em_px_emc", &_em_px_emc);
  _tree_KFP->Branch("_em_py_emc", &_em_py_emc);
  _tree_KFP->Branch("_em_pz_emc", &_em_pz_emc);
  _tree_KFP->Branch("_em_x_emc", &_em_x_emc);
  _tree_KFP->Branch("_em_y_emc", &_em_y_emc);
  _tree_KFP->Branch("_em_z_emc", &_em_z_emc);
  if(m_doSimulation)
  {
    _tree_KFP->Branch("_em_has_truthmatching", &_em_has_truthmatching);
    _tree_KFP->Branch("_em_true_id", &_em_true_id);
    _tree_KFP->Branch("_em_true_px", &_em_true_px);
    _tree_KFP->Branch("_em_true_py", &_em_true_py);
    _tree_KFP->Branch("_em_true_pz", &_em_true_pz);
    _tree_KFP->Branch("_em_true_vertex_x", &_em_true_vertex_x);
    _tree_KFP->Branch("_em_true_vertex_y", &_em_true_vertex_y);
    _tree_KFP->Branch("_em_true_vertex_z", &_em_true_vertex_z);
    _tree_KFP->Branch("_em_true_vertex_x_method2", &_em_true_vertex_x_method2);
    _tree_KFP->Branch("_em_true_vertex_y_method2", &_em_true_vertex_y_method2);
    _tree_KFP->Branch("_em_true_vertex_z_method2", &_em_true_vertex_z_method2);
  }

  _tree_KFP->Branch("_emcal_phi", &_emcal_phi);
  _tree_KFP->Branch("_emcal_eta", &_emcal_eta);
  _tree_KFP->Branch("_emcal_x", &_emcal_x);
  _tree_KFP->Branch("_emcal_y", &_emcal_y);
  _tree_KFP->Branch("_emcal_z", &_emcal_z);
  _tree_KFP->Branch("_emcal_e", &_emcal_e);

  _tree_KFP->Branch("_epem_DCA_2d", &_epem_DCA_2d);
  _tree_KFP->Branch("_epem_DCA_3d", &_epem_DCA_3d);

  if(m_doSimulation)
  {
    _tree_KFP->Branch("_true_numCan", &_true_numCan);

    _tree_KFP->Branch("_true_gamma_phi", &_true_gamma_phi);
    _tree_KFP->Branch("_true_gamma_eta", &_true_gamma_eta);
    _tree_KFP->Branch("_true_gamma_px", &_true_gamma_px);
    _tree_KFP->Branch("_true_gamma_py", &_true_gamma_py);
    _tree_KFP->Branch("_true_gamma_pz", &_true_gamma_pz);
    _tree_KFP->Branch("_true_gamma_pE", &_true_gamma_pE);
    _tree_KFP->Branch("_true_gamma_x", &_true_gamma_x);
    _tree_KFP->Branch("_true_gamma_y", &_true_gamma_y);
    _tree_KFP->Branch("_true_gamma_z", &_true_gamma_z);
    _tree_KFP->Branch("_true_gamma_mother_id", &_true_gamma_mother_id);
    _tree_KFP->Branch("_true_gamma_embedding_id", &_true_gamma_embedding_id);

    _tree_KFP->Branch("_true_ep_phi", &_true_ep_phi);
    _tree_KFP->Branch("_true_ep_eta", &_true_ep_eta);
    _tree_KFP->Branch("_true_ep_px", &_true_ep_px);
    _tree_KFP->Branch("_true_ep_py", &_true_ep_py);
    _tree_KFP->Branch("_true_ep_pz", &_true_ep_pz);
    _tree_KFP->Branch("_true_ep_pE", &_true_ep_pE);
    _tree_KFP->Branch("_true_ep_x", &_true_ep_x);
    _tree_KFP->Branch("_true_ep_y", &_true_ep_y);
    _tree_KFP->Branch("_true_ep_z", &_true_ep_z);

    _tree_KFP->Branch("_true_em_phi", &_true_em_phi);
    _tree_KFP->Branch("_true_em_eta", &_true_em_eta);
    _tree_KFP->Branch("_true_em_px", &_true_em_px);
    _tree_KFP->Branch("_true_em_py", &_true_em_py);
    _tree_KFP->Branch("_true_em_pz", &_true_em_pz);
    _tree_KFP->Branch("_true_em_pE", &_true_em_pE);
    _tree_KFP->Branch("_true_em_x", &_true_em_x);
    _tree_KFP->Branch("_true_em_y", &_true_em_y);
    _tree_KFP->Branch("_true_em_z", &_true_em_z);
   }
}

//____________________________________________________________________________..
int TrackToCalo::process_event(PHCompositeNode *topNode)
{
  PHNodeIterator nodeIter(topNode);
  PHNode* evtNode = dynamic_cast<PHNode*>(nodeIter.findFirst("EventHeader"));
  if (evtNode)
  {
    EventHeaderv1* evtHeader = findNode::getClass<EventHeaderv1>(topNode, "EventHeader");
    _runNumber = evtHeader->get_RunNumber();
    _eventNumber = evtHeader->get_EvtSequence();
  }
  else
  {
    _runNumber = 0;
    _eventNumber = -1;
  }
  std::cout<<"TrackToCalo::process_event run "<<_runNumber<<" event "<<_eventNumber<<std::endl;

  if(!trackMap)
  {
    trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
    if(!trackMap)
    {
      std::cout << "TrackToCalo::process_event: SvtxTrackMap not found!" << std::endl;
    }
  }

  if (!acts_Geometry)
  {
    acts_Geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
    if (!acts_Geometry)
    {
      std::cout << "TrackToCalo::process_event: ActsGeometry not found!" << std::endl;
    }
  }

  if (!clustersEM)
  {
    clustersEM = findNode::getClass<RawClusterContainer>(topNode, m_RawClusCont_EM_name);
    if (!clustersEM)
    {
      std::cout << "TrackToCalo::process_event: " << m_RawClusCont_EM_name << " not found!" << std::endl;
    }
  }
  if (!clustersHAD)
  {
    clustersHAD = findNode::getClass<RawClusterContainer>(topNode, m_RawClusCont_HAD_name);
    if (!clustersHAD)
    {
      std::cout << "TrackToCalo::process_event: " << m_RawClusCont_HAD_name << " not found!" << std::endl;
    }
  }

  if(!EMCAL_Container)
  {
    EMCAL_Container = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
    if(!EMCAL_Container)
    {
      std::cout << "TrackToCalo::process_event: TOWERINFO_CALIB_CEMC not found!" << std::endl;
    }
  }
  if(!IHCAL_Container)
  {
    IHCAL_Container = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
    if(!IHCAL_Container)
    {
      std::cout << "TrackToCalo::process_event: TOWERINFO_CALIB_HCALIN not found!" << std::endl;
    }
  }
  if(!OHCAL_Container)
  {
    OHCAL_Container = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
    if(!OHCAL_Container)
    {
      std::cout << "TrackToCalo::process_event: TOWERINFO_CALIB_HCALOUT not found!" << std::endl;
    }
  }

  if(!trkrContainer)
  {
    trkrContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    if(!trkrContainer)
    {
      std::cout << "TrackToCalo::process_event: TRKR_CLUSTER not found!" << std::endl;
    }
  }

  if(!EMCalGeo)
  {
    EMCalGeo = findNode::getClass<RawTowerGeomContainer>(topNode, m_RawTowerGeomCont_name);
    if(!EMCalGeo)
    {
      std::cout << "TrackToCalo::process_event: " << m_RawTowerGeomCont_name << " not found!" << std::endl;
    }
  }

  if(!IHCalGeo)
  {
    IHCalGeo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
    if(!IHCalGeo)
    {
      std::cout << "TrackToCalo::process_event: TOWERGEOM_HCALIN not found!" << std::endl;
    }
  }

  if(!OHCalGeo)
  {
    OHCalGeo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
    if(!OHCalGeo)
    {
      std::cout << "TrackToCalo::process_event: TOWERGEOM_HCALOUT not found!" << std::endl;
    }
  }

  if(!KFP_Container)
  {
    KFP_Container = findNode::getClass<KFParticle_Container>(topNode, m_KFPCont_name);
    if(!KFP_Container)
    {
      std::cout << "TrackToCalo::process_event: " << m_KFPCont_name << " not found!" << std::endl;
    }
  }

  if(!KFP_trackMap)
  {
    KFP_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_KFPtrackMap_name);
    if(!KFP_trackMap)
    {
      std::cout << "TrackToCalo::process_event: " << m_KFPtrackMap_name << " not found!" << std::endl;
    }
  }

  if(!vertexmap)
  {
    vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    if(!vertexmap)
    {
      std::cout << "TrackToCalo::process_event: GlobalVertexMap not found! (but not necessary)" << std::endl;
    }
  }

  if(!vertexMap)
  {
    vertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
    if(!vertexMap)
    {
      std::cout << "TrackToCalo::process_event: SvtxVertexMap not found! (but not necessary)" << std::endl;
    }
  }

  if(!gl1Packet)
  {
    gl1Packet = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    if(!gl1Packet)
    {
      std::cout << "TrackToCalo::process_event: GL1Packet not found! (but not necessary)" << std::endl;
    }
  }

  if(m_doSimulation)
  {
    if(!m_decayMap)
    {
      std::string df_node_name = m_df_module_name + "_DecayMap";
      m_decayMap = findNode::getClass<DecayFinderContainer_v1>(topNode, df_node_name.c_str());
      if(!m_decayMap)
      {
        std::cout << "TrackToCalo::process_event: " << df_node_name.c_str() << "not found! can not do truth matching" << std::endl;
      }
    }
  
    if (!m_geneventmap)
    {
      m_geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
      if (!m_geneventmap)
      {
        std::cout << "TrackToCalo::process_event: PHHepMCGenEventMap not found! can not do truth matching" << std::endl;
      }
    }
  
    if (!m_truthInfo)
    {
      m_truthInfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
      if (!m_truthInfo)
      {
        std::cout << "TrackToCalo::process_event: G4TruthInfo not found! can not do truth matching" << std::endl;
      }
    } 

    if(Verbosity() > 2)
    {
      std::cout << "Begin Print Truth Info" << std::endl;
      m_decayMap->identify();
      for (auto &iter : *m_decayMap)
      {
        Decay decay = iter.second;
  
        TLorentzVector *motherTrueLV = new TLorentzVector();
        TLorentzVector *daughterTrueLV = new TLorentzVector();
        TVector3 *mother3Vector = new TVector3();
        TVector3 *daughter3Vector = new TVector3();
  
        for (unsigned int i = 1; i < decay.size(); ++i)
        {
          int trackid = -999;
          int pid = -999;
          PHG4TruthInfoContainer::ConstRange range = m_truthInfo->GetParticleRange();
          for (PHG4TruthInfoContainer::ConstIterator iter2 = range.first; iter2 != range.second; ++iter2)
          {
            PHG4Particle *daughterG4 = iter2->second;
  
            PHG4Particle *motherG4 = nullptr;
            if (daughterG4->get_parent_id() != 0)
            {
              motherG4 = m_truthInfo->GetParticle(daughterG4->get_parent_id());
            }
            else
            {
              continue;
            }
  
            if (motherG4->get_pid() == decay[0].second && motherG4->get_barcode() == decay[0].first.second && daughterG4->get_pid() == decay[i].second && daughterG4->get_barcode() == decay[i].first.second)
            {
              pid = daughterG4->get_pid();
              trackid = daughterG4->get_track_id();
  
              TVector3 *motherTrue3Vector = new TVector3(motherG4->get_px(), motherG4->get_py(), motherG4->get_pz());
              motherTrueLV->SetVectM((*motherTrue3Vector), getParticleMass(decay[0].second));
  
              PHG4VtxPoint *thisVtx = m_truthInfo->GetVtx(motherG4->get_vtx_id());
              mother3Vector->SetXYZ(thisVtx->get_x(), thisVtx->get_y(), thisVtx->get_z());
  
              daughterTrueLV->SetVectM(TVector3(daughterG4->get_px(), daughterG4->get_py(), daughterG4->get_pz()), getParticleMass(decay[i].second));
  
              // Now get the decay vertex position
              thisVtx = m_truthInfo->GetVtx(daughterG4->get_vtx_id());
              daughter3Vector->SetXYZ(thisVtx->get_x(), thisVtx->get_y(), thisVtx->get_z());
  
              delete motherTrue3Vector;
            }
          }
  
  std::cout << " id " << pid
            << " trackid " << trackid
            << " px " << daughterTrueLV->Px()
            << " py " << daughterTrueLV->Py()
            << " pz " << daughterTrueLV->Pz()
            << " pt " << daughterTrueLV->Pt()
            << " p " << daughterTrueLV->P()
            << " e " << daughterTrueLV->E()
            << " eta " << daughterTrueLV->PseudoRapidity()
            << " phi " << daughterTrueLV->Phi()
            << " vx = " << daughter3Vector->X()
            << " vy = " << daughter3Vector->Y()
            << " vz = " << daughter3Vector->Z()
            << std::endl;
        }
  
  std::cout << " id " << 22
            << " px " << motherTrueLV->Px()
            << " py " << motherTrueLV->Py()
            << " pz " << motherTrueLV->Pz()
            << " pt " << motherTrueLV->Pt()
            << " p " << motherTrueLV->P()
            << " e " << motherTrueLV->E()
            << " eta " << motherTrueLV->PseudoRapidity()
            << " phi " << motherTrueLV->Phi()
            << " vx = " << mother3Vector->X()
            << " vy = " << mother3Vector->Y()
            << " vz = " << mother3Vector->Z()
            << std::endl;
      }
    }

    PHNode *svtxg4Node = dynamic_cast<PHNode *>(nodeIter.findFirst("SvtxPHG4ParticleMap"));
    //if (svtxg4Node && false)
    if (svtxg4Node)
    {
      svtxg4Node = dynamic_cast<PHNode *>(nodeIter.findFirst("G4TruthInfo"));
      if (svtxg4Node)
      {
        m_truthInfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
      }
      else
      {
        std::cout << "TrackToCalo truth matching: G4TruthInfo does not exist" << std::endl;
      }
      dst_reco_truth_map = findNode::getClass<SvtxPHG4ParticleMap_v1>(topNode, "SvtxPHG4ParticleMap");
    }
  
    if (m_truthInfo)
    {
      if (!m_svtx_evalstack)
      {
        m_svtx_evalstack = new SvtxEvalStack(topNode);
        // clustereval = m_svtx_evalstack->get_cluster_eval();
        // hiteval = m_svtx_evalstack->get_hit_eval();
        trackeval = m_svtx_evalstack->get_track_eval();
        trutheval = m_svtx_evalstack->get_truth_eval();
        vertexeval = m_svtx_evalstack->get_vertex_eval();
      }
      m_svtx_evalstack->next_event(topNode);
    }
  }

  if (m_doTrkrCaloMatching)
  {
    ResetTreeVectors();
    fillTree();
  }

  if (m_doTrkrCaloMatching_KFP)
  {
    ResetTreeVectors_KFP();
    fillTree_KFP();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void TrackToCalo::fillTree()
{
    if (m_doTrackOnly) {fillTree_TrackOnly();}
    if (m_doCaloOnly) {fillTree_CaloOnly();}
    if (m_doTrackOnly || m_doCaloOnly) {_tree->Fill();}
}

//____________________________________________________________________________..
void TrackToCalo::fillTree_TrackOnly()
{
  if (!trackMap || !acts_Geometry || !trkrContainer)
  {
    std::cout << PHWHERE << "missing node trees, can't continue with track-calo matching (Track Only Part)"
              << std::endl;
    return;
  }

  bool has_vertex = false;
  GlobalVertex *mbd_vtx = nullptr;

  CLHEP::Hep3Vector vertex(0., 0., 0.);

  if(vertexmap)
  {
    if(!vertexmap->empty())
    {
      mbd_vtx = vertexmap->begin()->second;
      if(mbd_vtx)
      {
        vertex.setX(mbd_vtx->get_x());
        vertex.setY(mbd_vtx->get_y());
        vertex.setZ(mbd_vtx->get_z());
        has_vertex = true;
      }
    }
  }

  if(has_vertex)
  {
    _mbd_x.push_back(vertex.x());
    _mbd_y.push_back(vertex.y());
    _mbd_z.push_back(vertex.z());
  }
  else
  {
    _mbd_x.push_back(NAN);
    _mbd_y.push_back(NAN);
    _mbd_z.push_back(NAN);
  }

  if(vertexMap)
  {
    if(!vertexMap->empty())
    {
      for (const auto& [key, svtx_vtx] : *vertexMap)
      {
        _vertex_id.push_back(svtx_vtx->get_id());
        _vertex_crossing.push_back(svtx_vtx->get_beam_crossing());
        _vertex_ntracks.push_back(svtx_vtx->size_tracks());
        _vertex_x.push_back(svtx_vtx->get_x());
        _vertex_y.push_back(svtx_vtx->get_y());
        _vertex_z.push_back(svtx_vtx->get_z());
      }
    }
  }

  if(gl1Packet)
  {
    auto scaled_vector = gl1Packet->getScaledVector();
    for(int i = 0; i < 32; i++)
    {
      if((scaled_vector & (int)std::pow(2,i)) != 0)
      {
        _triggers.push_back(i);
      }
    }
  }

  _ntracks.push_back(trackMap->size());

  TrkrClusterContainer::HitSetKeyList tpcHits = trkrContainer->getHitSetKeys(TrkrDefs::TrkrId::tpcId);
  for (auto &hsk : tpcHits)
  {
    auto range = trkrContainer->getClusters(hsk);
    for (auto iter = range.first; iter != range.second; ++iter)
    {
      const auto cluskey = iter->first;
      const auto cluster = iter->second;  // auto cluster = clusters->findCluster(key);
      auto glob = acts_Geometry->getGlobalPosition(cluskey, cluster);
      auto sclusgx = glob.x();
      auto sclusgy = glob.y();
      auto sclusgz = glob.z();
      _cluster_x.push_back(sclusgx);
      _cluster_y.push_back(sclusgy);
      _cluster_z.push_back(sclusgz);
    }
  }

  for (auto &iter : *trackMap)
  {
    svtxtrack = iter.second;

    if(!checkTrack(svtxtrack))
    {
      continue;
    }

    si_seed = svtxtrack->get_silicon_seed();
    tpc_seed = svtxtrack->get_tpc_seed();

    int n_mvtx_clusters = 0;
    int n_intt_clusters = 0;
    int n_tpc_clusters = 0;
    int n_tpot_clusters = 0;

    if(!si_seed)
    {
      _track_nc_mvtx.push_back(0);
      _track_nc_intt.push_back(0);
      _track_bc.push_back(SHRT_MAX);
    }
    else
    {
      for(auto key_iter = si_seed->begin_cluster_keys(); key_iter != si_seed->end_cluster_keys(); ++key_iter)
      {
        const auto& cluster_key = *key_iter;
        trkrCluster = trkrContainer->findCluster(cluster_key);
        if(!trkrCluster)
        {
          continue;
        }
        if(TrkrDefs::getTrkrId(cluster_key) == TrkrDefs::TrkrId::mvtxId)
        {
          n_mvtx_clusters++;
        }
        if(TrkrDefs::getTrkrId(cluster_key) == TrkrDefs::TrkrId::inttId)
        {
          n_intt_clusters++;
        }
        Acts::Vector3 global(0., 0., 0.);
        global = acts_Geometry->getGlobalPosition(cluster_key, trkrCluster);
        _trClus_track_id.push_back(svtxtrack->get_id());
        _trClus_type.push_back(TrkrDefs::getTrkrId(cluster_key));
        _trClus_x.push_back(global[0]);
        _trClus_y.push_back(global[1]);
        _trClus_z.push_back(global[2]);
      }
      _track_nc_mvtx.push_back(n_mvtx_clusters);
      _track_nc_intt.push_back(n_intt_clusters);
      _track_bc.push_back(si_seed->get_crossing());
    }

    if(tpc_seed)
    {
      for(auto key_iter = tpc_seed->begin_cluster_keys(); key_iter != tpc_seed->end_cluster_keys(); ++key_iter)
      {
        const auto& cluster_key = *key_iter;
        trkrCluster = trkrContainer->findCluster(cluster_key);
        if(!trkrCluster)
        {
          continue;
        }
        if(TrkrDefs::getTrkrId(cluster_key) == TrkrDefs::TrkrId::tpcId)
        {
          n_tpc_clusters++;
        }
        if(TrkrDefs::getTrkrId(cluster_key) == TrkrDefs::TrkrId::micromegasId)
        {
          n_tpot_clusters++;
        }
        Acts::Vector3 global(0., 0., 0.);
        global = acts_Geometry->getGlobalPosition(cluster_key, trkrCluster);
        _trClus_track_id.push_back(svtxtrack->get_id());
        _trClus_type.push_back(TrkrDefs::getTrkrId(cluster_key));
        _trClus_x.push_back(global[0]);
        _trClus_y.push_back(global[1]);
        _trClus_z.push_back(global[2]);
      }
      _track_nc_tpc.push_back(n_tpc_clusters);
      _track_nc_tpot.push_back(n_tpot_clusters);
    }

    unsigned int m_vertexid = svtxtrack->get_vertex_id();
    bool track_have_vertex = false;
    if (vertexMap)
    {
      auto vertexit = vertexMap->find(m_vertexid);
      if (vertexit != vertexMap->end())
      {
        auto svtxvertex = vertexit->second;
        m_vx = svtxvertex->get_x();
        m_vy = svtxvertex->get_y();
        m_vz = svtxvertex->get_z();
        track_have_vertex = true;
      }
    }

    if (track_have_vertex)
    {
      _track_vx.push_back(m_vx);
      _track_vy.push_back(m_vy);
      _track_vz.push_back(m_vz);
    }
    else
    {
      _track_vx.push_back(NAN);
      _track_vy.push_back(NAN);
      _track_vz.push_back(NAN);
    }

    _track_id.push_back(svtxtrack->get_id());
    _track_quality.push_back(svtxtrack->get_quality());
    //auto dcapair = TrackAnalysisUtils::get_dca(svtxtrack, acts_vertex);
    Acts::Vector3 zero = Acts::Vector3::Zero();
    auto dcapair = TrackAnalysisUtils::get_dca(svtxtrack, zero);
    _track_dcaxy.push_back(dcapair.first.first);
    _track_dcaz.push_back(dcapair.second.first);
    _track_ptq.push_back(svtxtrack->get_charge()*svtxtrack->get_pt());
    _track_px.push_back(svtxtrack->get_px());
    _track_py.push_back(svtxtrack->get_py());
    _track_pz.push_back(svtxtrack->get_pz());
    _track_phi.push_back(svtxtrack->get_phi());
    _track_eta.push_back(svtxtrack->get_eta());
    _track_pcax.push_back(svtxtrack->get_x());
    _track_pcay.push_back(svtxtrack->get_y());
    _track_pcaz.push_back(svtxtrack->get_z());
    _track_crossing.push_back(svtxtrack->get_crossing());

    if (!m_doCaloOnly) continue;

    resetCaloRadius();

    // project to R=0
    thisState = svtxtrack->get_state(0);

    if(!thisState)
    {
      _track_phi_origin.push_back(NAN);
      _track_eta_origin.push_back(NAN);
      _track_px_origin.push_back(NAN);
      _track_py_origin.push_back(NAN);
      _track_pz_origin.push_back(NAN);
      _track_x_origin.push_back(NAN);
      _track_y_origin.push_back(NAN);
      _track_z_origin.push_back(NAN);
    }
    else
    {
      _track_phi_origin.push_back(thisState->get_phi());
      _track_eta_origin.push_back(thisState->get_eta());
      _track_px_origin.push_back(thisState->get_px());
      _track_py_origin.push_back(thisState->get_py());
      _track_pz_origin.push_back(thisState->get_pz());
      _track_x_origin.push_back(thisState->get_x());
      _track_y_origin.push_back(thisState->get_y());
      _track_z_origin.push_back(thisState->get_z());
    }

    // project to R_EMCAL
    thisState = svtxtrack->get_state(caloRadiusEMCal);

    if(!thisState)
    {
      _track_phi_emc.push_back(NAN);
      _track_eta_emc.push_back(NAN);
      _track_px_emc.push_back(NAN);
      _track_py_emc.push_back(NAN);
      _track_pz_emc.push_back(NAN);
      _track_x_emc.push_back(NAN);
      _track_y_emc.push_back(NAN);
      _track_z_emc.push_back(NAN);
    }
    else
    {
      _track_phi_emc.push_back(thisState->get_phi());
      _track_eta_emc.push_back(thisState->get_eta());
      _track_px_emc.push_back(thisState->get_px());
      _track_py_emc.push_back(thisState->get_py());
      _track_pz_emc.push_back(thisState->get_pz());
      _track_x_emc.push_back(thisState->get_x());
      _track_y_emc.push_back(thisState->get_y());
      _track_z_emc.push_back(thisState->get_z());
    }

    // project to R_IHCAL
    thisState = svtxtrack->get_state(caloRadiusIHCal);

    if(!thisState)
    {
      _track_phi_ihc.push_back(NAN);
      _track_eta_ihc.push_back(NAN);
      _track_px_ihc.push_back(NAN);
      _track_py_ihc.push_back(NAN);
      _track_pz_ihc.push_back(NAN);
      _track_x_ihc.push_back(NAN);
      _track_y_ihc.push_back(NAN);
      _track_z_ihc.push_back(NAN);
    }
    else
    {
      _track_phi_ihc.push_back(thisState->get_phi());
      _track_eta_ihc.push_back(thisState->get_eta());
      _track_px_ihc.push_back(thisState->get_px());
      _track_py_ihc.push_back(thisState->get_py());
      _track_pz_ihc.push_back(thisState->get_pz());
      _track_x_ihc.push_back(thisState->get_x());
      _track_y_ihc.push_back(thisState->get_y());
      _track_z_ihc.push_back(thisState->get_z());
    }

    // project to R_OHCAL
    thisState = svtxtrack->get_state(caloRadiusOHCal);

    if(!thisState)
    {
      _track_phi_ohc.push_back(NAN);
      _track_eta_ohc.push_back(NAN);
      _track_px_ohc.push_back(NAN);
      _track_py_ohc.push_back(NAN);
      _track_pz_ohc.push_back(NAN);
      _track_x_ohc.push_back(NAN);
      _track_y_ohc.push_back(NAN);
      _track_z_ohc.push_back(NAN);
    }
    else
    {
      _track_phi_ohc.push_back(thisState->get_phi());
      _track_eta_ohc.push_back(thisState->get_eta());
      _track_px_ohc.push_back(thisState->get_px());
      _track_py_ohc.push_back(thisState->get_py());
      _track_pz_ohc.push_back(thisState->get_pz());
      _track_x_ohc.push_back(thisState->get_x());
      _track_y_ohc.push_back(thisState->get_y());
      _track_z_ohc.push_back(thisState->get_z());
    }
  }
}

//____________________________________________________________________________..
void TrackToCalo::fillTree_CaloOnly()
{
  if (!clustersEM || !clustersHAD || !EMCAL_Container || !IHCAL_Container || !OHCAL_Container || !EMCalGeo || !IHCalGeo || !OHCalGeo)
  {
    std::cout << PHWHERE << "missing node trees, can't continue with track-calo matching (Calo Only Part)"
              << std::endl;
    return;
  }

  bool has_vertex = false;
  GlobalVertex *mbd_vtx = nullptr;

  CLHEP::Hep3Vector vertex(0., 0., 0.);

  if(vertexmap)
  {
    if(!vertexmap->empty())
    {
      mbd_vtx = vertexmap->begin()->second;
      if(mbd_vtx)
      {
        vertex.setX(mbd_vtx->get_x());
        vertex.setY(mbd_vtx->get_y());
        vertex.setZ(mbd_vtx->get_z());
        has_vertex = true;
      }
    }
  }

  if(!m_doTrackOnly)
  {
    if(has_vertex)
    {
      _mbd_x.push_back(vertex.x());
      _mbd_y.push_back(vertex.y());
      _mbd_z.push_back(vertex.z());
    }
    else
    {
      _mbd_x.push_back(NAN);
      _mbd_y.push_back(NAN);
      _mbd_z.push_back(NAN);
    }
  }

  if(gl1Packet && !m_doTrackOnly)
  {
    auto scaled_vector = gl1Packet->getScaledVector();
    for(int i = 0; i < 32; i++)
    {
      if((scaled_vector & (int)std::pow(2,i)) != 0)
      {
        _triggers.push_back(i);
      }
    }
  }

  /*
  resetCaloRadius();

  caloRadiusEMCal *= Acts::UnitConstants::cm;

  const auto eta = 2.5;
  const auto theta = 2. * atan(exp(-eta));
  const auto halfZ = caloRadiusEMCal / tan(theta) * Acts::UnitConstants::cm;

  auto transform = Acts::Transform3::Identity();

  std::shared_ptr<Acts::CylinderSurface> emcal_surf = Acts::Surface::makeShared<Acts::CylinderSurface>(transform, caloRadiusEMCal, halfZ);


  caloRadiusIHCal *= Acts::UnitConstants::cm;

  const auto halfZ_IHCal = caloRadiusIHCal / tan(theta) * Acts::UnitConstants::cm;

  auto transform_IHCal = Acts::Transform3::Identity();

  std::shared_ptr<Acts::CylinderSurface> ihcal_surf = Acts::Surface::makeShared<Acts::CylinderSurface>(transform_IHCal, caloRadiusIHCal, halfZ_IHCal);

  caloRadiusOHCal *= Acts::UnitConstants::cm;

  const auto halfZ_OHCal = caloRadiusOHCal / tan(theta) * Acts::UnitConstants::cm;

  auto transform_OHCal = Acts::Transform3::Identity();

  std::shared_ptr<Acts::CylinderSurface> ohcal_surf = Acts::Surface::makeShared<Acts::CylinderSurface>(transform_OHCal, caloRadiusOHCal, halfZ_OHCal);


  TrackSeed *seed = nullptr;

  RawCluster *cluster = nullptr;

  ActsPropagator prop(acts_Geometry);

  TowerInfo *tInfo = nullptr;

  for(std::size_t i = 0; i < seedContainer->size(); i++)
  {
    seed = seedContainer->get(i);
    if(!seed)
    {
      continue;
    }

    if(seed->get_pt() < 0.5)
    {
      continue;
    }

    TrkrCluster *trkrCluster = nullptr;

    Acts::Vector3 global(0., 0., 0.);

    TrkrDefs::cluskey c_key = UINT64_MAX;

    for (auto key_iter = seed->begin_cluster_keys(); key_iter != seed->end_cluster_keys(); ++key_iter)
    {
      const auto& cluster_key = *key_iter;

      trkrCluster = trkrContainer->findCluster(cluster_key);

      c_key = cluster_key;

      global = acts_Geometry->getGlobalPosition(cluster_key, trkrCluster);

    }

    std::shared_ptr<const Acts::Surface> surf_mvtx = acts_Geometry->maps().getSurface(c_key, trkrCluster);


    SvtxTrackState_v1 *state = new SvtxTrackState_v1(0.);
    state->set_x(global[0]);
    state->set_y(global[1]);
    state->set_z(global[2]);
    float c_phi = std::atan2(global[1], global[0]);
    float c_r = sqrt(global[0]*global[0] + global[1]*global[1]);
    float c_theta = std::atan2(c_r, global[2]);
    float c_eta = -std::log(std::tan(c_theta / 2.));

    float c_px = seed->get_pt()*std::cos(c_phi);
    float c_py = seed->get_pt()*std::sin(c_phi);
    float c_pz = seed->get_pt()*std::sinh(c_eta);

    state->set_px(c_px);
    state->set_py(c_py);
    state->set_pz(c_pz);

    auto params = prop.makeTrackParams(state, seed->get_charge(), surf_mvtx);

    if(!params.ok())
    {
      std::cout << "NOT OK!" << std::endl;
      continue;
    }

    ActsPropagator propagator(acts_Geometry);
    propagator.constField();
    propagator.verbosity(Verbosity());
    propagator.setConstFieldValue(1.4 * Acts::UnitConstants::T);

    auto result_EMCal = propagator.propagateTrackFast((const Acts::BoundTrackParameters&)params, emcal_surf);

    auto result_IHCal = propagator.propagateTrackFast((const Acts::BoundTrackParameters&)params, ihcal_surf);

    auto result_OHCal = propagator.propagateTrackFast((const Acts::BoundTrackParameters&)params, ohcal_surf);

    float EMCal_proj_phi = 9999.;
    float IHCal_proj_phi = 9999.;
    float OHCal_proj_phi = 9999.;

    if (result_EMCal.ok())
    {
      auto projectionPos = result_EMCal.value().second.position(acts_Geometry->geometry().getGeoContext());
      EMCal_proj_phi = std::atan2(projectionPos.y(),projectionPos.x());
    }
    else
    {
      std::cout << "result_EMCal failed!" << std::endl;
      continue;
    }

    if (result_IHCal.ok())
    {
      auto projectionPos = result_IHCal.value().second.position(acts_Geometry->geometry().getGeoContext());
      IHCal_proj_phi = std::atan2(projectionPos.y(),projectionPos.x());
    }
    else
    {
      std::cout << "result_IHCal failed!" << std::endl;
      continue;
    }

    if (result_OHCal.ok())
    {
      auto projectionPos = result_OHCal.value().second.position(acts_Geometry->geometry().getGeoContext());
      OHCal_proj_phi = std::atan2(projectionPos.y(),projectionPos.x());
    }
    else
    {
      std::cout << "result_OHCal failed!" << std::endl;
      continue;
    }

    if(seed->get_charge() > 0) _EMCalProjShiftPhiPos->Fill(seed->get_phi()-EMCal_proj_phi);
    if(seed->get_charge() < 0) _EMCalProjShiftPhiNeg->Fill(seed->get_phi()-EMCal_proj_phi);

    _EMCalProjShiftPhiVsSeedPt->Fill(seed->get_pt(), seed->get_phi()-EMCal_proj_phi);

    _SeedPt->Fill(seed->get_pt());



    track_mult += 1.;

    _SeedPhi->Fill(seed->get_phi());

    bool EMCal_match = false;
    bool IHCal_match = false;
    bool OHCal_match = false;

    RawClusterContainer::Range begin_end_EMC = EMCAL_RawClusters->getClusters();
    RawClusterContainer::Iterator clusIter_EMC;

    /// Loop over the EMCal clusters
    for (clusIter_EMC = begin_end_EMC.first; clusIter_EMC != begin_end_EMC.second; ++clusIter_EMC)
    {
      cluster = clusIter_EMC->second;
      if(cluster->get_energy() < 0.2) continue;



      float deltaPhi = EMCal_proj_phi - RawClusterUtility::GetAzimuthAngle(*cluster, vertex);
      if(deltaPhi > M_PI)
      {
        deltaPhi = deltaPhi - 2*M_PI;
      }
      if(deltaPhi < -M_PI)
      {
        deltaPhi = deltaPhi + 2*M_PI;
      }

      _TrackEMCalPhi->Fill(deltaPhi);

      if(std::fabs(deltaPhi) < 0.2)
      {
        _TrackEMCalMatchPhi->Fill(deltaPhi);
        EMCal_match = true;
      }
    }

    //IHCal Towers
    for(unsigned int ihcal = 0; ihcal < IHCAL_Container->size(); ihcal++)
    {
      tInfo = IHCAL_Container->get_tower_at_channel(ihcal);

      if(!tInfo)
      {
        continue;
      }

      if(tInfo->get_energy() < 0.2) continue;

      unsigned int towerinfo_key = IHCAL_Container->encode_key(ihcal);
      int ti_ieta = IHCAL_Container->getTowerEtaBin(towerinfo_key);
      int ti_iphi = IHCAL_Container->getTowerPhiBin(towerinfo_key);
      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ti_ieta, ti_iphi);

      RawTowerGeom *tower_geom = IHCalGeo->get_tower_geometry(key);

      float deltaPhi = IHCal_proj_phi - tower_geom->get_phi();

      if(deltaPhi > M_PI)
      {
        deltaPhi = deltaPhi - 2*M_PI;
      }
      if(deltaPhi < -M_PI)
      {
        deltaPhi = deltaPhi + 2*M_PI;
      }

      _TrackIHCalPhi->Fill(deltaPhi);

      if(fabs(deltaPhi) < 0.2)
      {
        IHCal_match = true;
      }

    }

    //OHCal Towers
    for(unsigned int ohcal = 0; ohcal < OHCAL_Container->size(); ohcal++)
    {
      tInfo = OHCAL_Container->get_tower_at_channel(ohcal);

      if(!tInfo)
      {
        continue;
      }

      if(tInfo->get_energy() < 0.2) continue;

      unsigned int towerinfo_key = OHCAL_Container->encode_key(ohcal);
      int ti_ieta = OHCAL_Container->getTowerEtaBin(towerinfo_key);
      int ti_iphi = OHCAL_Container->getTowerPhiBin(towerinfo_key);
      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, ti_ieta, ti_iphi);

      RawTowerGeom *tower_geom = OHCalGeo->get_tower_geometry(key);

      float deltaPhi = OHCal_proj_phi - tower_geom->get_phi();

      if(deltaPhi > M_PI)
      {
        deltaPhi = deltaPhi - 2*M_PI;
      }
      if(deltaPhi < -M_PI)
      {
        deltaPhi = deltaPhi + 2*M_PI;
      }

      _TrackOHCalPhi->Fill(deltaPhi);

      if(fabs(deltaPhi) < 0.2)
      {
        OHCal_match = true;
      }


    }

  }
  */

  RawCluster *cluster = nullptr;

  RawClusterContainer::Range begin_end_EMC = clustersEM->getClusters();
  RawClusterContainer::Iterator clusIter_EMC;

  /// Loop over the EMCal clusters
  for (clusIter_EMC = begin_end_EMC.first; clusIter_EMC != begin_end_EMC.second; ++clusIter_EMC)
  {
    cluster = clusIter_EMC->second;
    if(cluster->get_energy() < m_emcal_e_low_cut) continue;

    _emcal_id.push_back(clusIter_EMC->first);
    _emcal_e.push_back(cluster->get_energy());
    _emcal_phi.push_back(RawClusterUtility::GetAzimuthAngle(*cluster, vertex));
    _emcal_eta.push_back(RawClusterUtility::GetPseudorapidity(*cluster, vertex));
    _emcal_x.push_back(cluster->get_x());
    _emcal_y.push_back(cluster->get_y());
    _emcal_z.push_back(cluster->get_z());
    _emcal_ecore.push_back(cluster->get_ecore());
    _emcal_chi2.push_back(cluster->get_chi2());
    _emcal_prob.push_back(cluster->get_prob());

    RawCluster::TowerConstRange towers = cluster->get_towers();
    RawCluster::TowerConstIterator toweriter;

    TowerInfo *towerInfo = nullptr;

    for (toweriter = towers.first; toweriter != towers.second; ++toweriter)
    {
      _emcal_tower_cluster_id.push_back(clusIter_EMC->first);
      RawTowerGeom *tower_geom = EMCalGeo->get_tower_geometry(toweriter->first);
      _emcal_tower_phi.push_back(tower_geom->get_phi());
      _emcal_tower_eta.push_back(tower_geom->get_eta());
      _emcal_tower_e.push_back(toweriter->second);
      unsigned int key = TowerInfoDefs::encode_emcal(tower_geom->get_bineta(), tower_geom->get_binphi());
      towerInfo = EMCAL_Container->get_tower_at_key(key);
      _emcal_tower_status.push_back(towerInfo->get_status());

    }
  }

  RawClusterContainer::Range begin_end_HAD = clustersHAD->getClusters();
  RawClusterContainer::Iterator clusIter_HAD;

  /// Loop over the HCal clusters
  for (clusIter_HAD = begin_end_HAD.first; clusIter_HAD != begin_end_HAD.second; ++clusIter_HAD)
  {
    cluster = clusIter_HAD->second;
    if(cluster->get_energy() < m_hcal_e_low_cut) continue;

    _hcal_id.push_back(clusIter_HAD->first);
    _hcal_e.push_back(cluster->get_energy());
    _hcal_phi.push_back(RawClusterUtility::GetAzimuthAngle(*cluster, vertex));
    _hcal_eta.push_back(RawClusterUtility::GetPseudorapidity(*cluster, vertex));
    _hcal_x.push_back(cluster->get_x());
    _hcal_y.push_back(cluster->get_y());
    _hcal_z.push_back(cluster->get_z());

    RawCluster::TowerConstRange towers = cluster->get_towers();
    RawCluster::TowerConstIterator toweriter;

    for (toweriter = towers.first; toweriter != towers.second; ++toweriter)
    {
      _hcal_tower_cluster_id.push_back(clusIter_HAD->first);
      RawTowerGeom *tower_geom = nullptr;
      TowerInfo *towerInfo = nullptr;
      int tower_io = -1;

      if(RawTowerDefs::decode_caloid(toweriter->first) == RawTowerDefs::CalorimeterId::HCALOUT)
      {
        tower_geom = OHCalGeo->get_tower_geometry(toweriter->first);
        unsigned int key = TowerInfoDefs::encode_hcal(tower_geom->get_bineta(), tower_geom->get_binphi());
        towerInfo = OHCAL_Container->get_tower_at_key(key);
        tower_io = 2;
      }
      else if(RawTowerDefs::decode_caloid(toweriter->first) == RawTowerDefs::CalorimeterId::HCALIN)
      {
        tower_geom = IHCalGeo->get_tower_geometry(toweriter->first);
        unsigned int key = TowerInfoDefs::encode_hcal(tower_geom->get_bineta(), tower_geom->get_binphi());
        towerInfo = IHCAL_Container->get_tower_at_key(key);
        tower_io = 1;
      }

      _hcal_tower_phi.push_back(tower_geom->get_phi());
      _hcal_tower_eta.push_back(tower_geom->get_eta());
      _hcal_tower_e.push_back(toweriter->second);
      _hcal_tower_status.push_back(towerInfo->get_status());
      _hcal_tower_io.push_back(tower_io);
    }
  }

  /*

  RawClusterContainer::Range begin_end_EMC = EMCAL_RawClusters->getClusters();
  RawClusterContainer::Iterator clusIter_EMC;

  /// Loop over the EMCal clusters
  for (clusIter_EMC = begin_end_EMC.first; clusIter_EMC != begin_end_EMC.second; ++clusIter_EMC)
  {
    cluster = clusIter_EMC->second;
    if(cluster->get_energy() < 0.2) continue;

    _emcal_e.push_back(cluster->get_energy());
    _emcal_phi.push_back(RawClusterUtility::GetAzimuthAngle(*cluster, vertex));
    _emcal_eta.push_back(RawClusterUtility::GetPseudorapidity(*cluster, vertex));

  }

  for(unsigned int ihcal = 0; ihcal < IHCAL_Container->size(); ihcal++)
  {
    tInfo = IHCAL_Container->get_tower_at_channel(ihcal);

    if(!tInfo)
    {
      continue;
    }

    if(tInfo->get_energy() < 0.2) continue;

    unsigned int towerinfo_key = IHCAL_Container->encode_key(ihcal);
    int ti_ieta = IHCAL_Container->getTowerEtaBin(towerinfo_key);
    int ti_iphi = IHCAL_Container->getTowerPhiBin(towerinfo_key);
    const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ti_ieta, ti_iphi);

    RawTowerGeom *tower_geom = IHCalGeo->get_tower_geometry(key);

    _ihcal_e.push_back(tInfo->get_energy());
    _ihcal_phi.push_back(tower_geom->get_phi());
    _ihcal_eta.push_back(tower_geom->get_eta());

  }

  for(unsigned int ohcal = 0; ohcal < OHCAL_Container->size(); ohcal++)
  {
    tInfo = OHCAL_Container->get_tower_at_channel(ohcal);

    if(!tInfo)
    {
      continue;
    }

    if(tInfo->get_energy() < 0.2) continue;

    unsigned int towerinfo_key = OHCAL_Container->encode_key(ohcal);
    int ti_ieta = OHCAL_Container->getTowerEtaBin(towerinfo_key);
    int ti_iphi = OHCAL_Container->getTowerPhiBin(towerinfo_key);
    const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, ti_ieta, ti_iphi);

    RawTowerGeom *tower_geom = OHCalGeo->get_tower_geometry(key);

    _ohcal_e.push_back(tInfo->get_energy());
    _ohcal_phi.push_back(tower_geom->get_phi());
    _ohcal_eta.push_back(tower_geom->get_eta());

  }

  */
}

void TrackToCalo::fillTree_KFP()
{
  if (!KFP_Container || !KFP_trackMap || !acts_Geometry || !clustersEM || !EMCalGeo || !trkrContainer)
  {
    std::cout << PHWHERE << "missing node trees, can't continue with track-calo matching with KFParticle"
              << std::endl;
    return;
  }

  CLHEP::Hep3Vector vertex(0., 0., 0.);

  resetCaloRadius();

  if (KFP_Container->empty())
  {
    std::cout<<"No KFParticle reconstructed in this event! Skip!"<<std::endl;
    return;
  }

  size_t length_kfps = KFP_Container->size();
  if (static_cast<int>(length_kfps) % 3 != 0)
  {
    std::cout<<"Why KFParticle is not 3*n? Skip!"<<std::endl;
    return;
  }

  if (Verbosity() > 2)
  {
    for (auto &iter : *KFP_Container)
    {
      KFParticle* kfp = iter.second;
      std::cout<<"KFP ID = "<<kfp->Id()<<" PDGID = "<<kfp->GetPDG()<<" p = "<<kfp->GetP()<<std::endl;
    }

    for (auto &iter : *KFP_trackMap)
    {
      SvtxTrack* kfp = iter.second;
      std::cout<<"iter.first = "<<iter.first<<std::endl;
      std::cout<<"svtxtrack id = "<<kfp->get_id()<<" px = "<<kfp->get_px()<<" py = "<<kfp->get_py()<<" pz = "<<kfp->get_pz()<<std::endl;
    }
  }

  _numCan = static_cast<int>(length_kfps) / 3;

  for (int i = 0; i < _numCan; i++)
  {
    auto it_kfp_cont = KFP_Container->begin();
    std::advance(it_kfp_cont, 3 * i);
    kfp_mother = it_kfp_cont->second;

    float mass, massErr;
    kfp_mother->GetMass(mass, massErr);
    _gamma_mass.push_back(mass);
    _gamma_massErr.push_back(massErr);
    _gamma_x.push_back(kfp_mother->GetX());
    _gamma_y.push_back(kfp_mother->GetY());
    _gamma_z.push_back(kfp_mother->GetZ());
    _gamma_px.push_back(kfp_mother->GetPx());
    _gamma_py.push_back(kfp_mother->GetPy());
    _gamma_pz.push_back(kfp_mother->GetPz());
    _gamma_pE.push_back(kfp_mother->GetE());
    _gamma_pT.push_back(kfp_mother->GetPt());
    _gamma_pTErr.push_back(kfp_mother->GetErrPt());
    _gamma_p.push_back(kfp_mother->GetP());
    _gamma_pErr.push_back(kfp_mother->GetErrP());
    _gamma_pseudorapidity.push_back(kfp_mother->GetEta());
    _gamma_rapidity.push_back(kfp_mother->GetRapidity());
    _gamma_theta.push_back(kfp_mother->GetTheta());
    _gamma_phi.push_back(kfp_mother->GetPhi());
    _gamma_chi2.push_back(kfp_mother->GetChi2());
    _gamma_nDoF.push_back(kfp_mother->GetNDF());
    _gamma_vertex_volume.push_back( kf_tools.calculateEllipsoidVolume(*kfp_mother) );

    int totalcharge = 0;
    // one for e+, one for e-
    for (int j = 1; j <= 2; j++)
    {

      it_kfp_cont = KFP_Container->begin();
      std::advance(it_kfp_cont, 3 * i + j);
      kfp_daughter = it_kfp_cont->second;
      float p_daughter_unmoved = kfp_daughter->GetP();
      float e_daughter_unmoved = kfp_daughter->GetE();
      float pt_daughter_unmoved = kfp_daughter->GetPt();

      if (Verbosity() > 1)
      {
        std::cout<<"KFP: before SetProductionVertex kfp_daughter->GetMass() = "<<kfp_daughter->GetMass()<<std::endl;
        std::cout<<"KFP: before SetProductionVertex sqrt(E^2-p^2) = "<<sqrt(pow(kfp_daughter->GetE(),2) - pow(kfp_daughter->GetPx(),2) - pow(kfp_daughter->GetPy(),2) - pow(kfp_daughter->GetPz(),2))<<std::endl;
      }

      kfp_daughter->SetProductionVertex(*kfp_mother);

      if (Verbosity() > 1)
      {
        std::cout<<"KFP: after SetProductionVertex kfp_daughter->GetMass() = "<<kfp_daughter->GetMass()<<std::endl;
        std::cout<<"KFP: after SetProductionVertex sqrt(E^2-p^2) = "<<sqrt(pow(kfp_daughter->GetE(),2) - pow(kfp_daughter->GetPx(),2) - pow(kfp_daughter->GetPy(),2) - pow(kfp_daughter->GetPz(),2))<<std::endl;
      }

      auto it_kfp_trackmap = KFP_trackMap->begin();
      std::advance(it_kfp_trackmap, 3 * i + j);
      svtxtrack = it_kfp_trackmap->second;

      bool isParticleValid = false;
      int true_id = 0;
      float true_px = 0, true_py = 0, true_pz = 0;
      float true_daughter_vertex_x = 0, true_daughter_vertex_y = 0, true_daughter_vertex_z = 0;
      float true_daughter_vertex_x_method2 = 0, true_daughter_vertex_y_method2 = 0, true_daughter_vertex_z_method2 = 0;
      if (m_doSimulation && m_doTruthMatching)
      {
        g4particle = getTruthTrack(svtxtrack);

        isParticleValid = g4particle == nullptr ? false : true;

        true_id = isParticleValid ? g4particle->get_pid() : 0;
        true_px = isParticleValid ? (Float_t) g4particle->get_px() : 0.;
        true_py = isParticleValid ? (Float_t) g4particle->get_py() : 0.;
        true_pz = isParticleValid ? (Float_t) g4particle->get_pz() : 0.;

        if (isParticleValid)
        {
          g4vertex_point = trutheval->get_vertex(g4particle);
        }

        true_daughter_vertex_x = isParticleValid ? g4vertex_point->get_x() : 0.;
        true_daughter_vertex_y = isParticleValid ? g4vertex_point->get_y() : 0.;
        true_daughter_vertex_z = isParticleValid ? g4vertex_point->get_z() : 0.;

        if (g4particle)
        {
          // Now get the decay vertex position
          PHG4VtxPoint *thisVtx = m_truthInfo->GetVtx(g4particle->get_vtx_id());
          true_daughter_vertex_x_method2 = thisVtx->get_x();
          true_daughter_vertex_y_method2 = thisVtx->get_y();
          true_daughter_vertex_z_method2 = thisVtx->get_z();
	  if (Verbosity() > 2)
	  {
            std::cout<<"g4particle->get_pid() = "<<g4particle->get_pid()<<" , g4particle->get_barcode() = "<<g4particle->get_barcode()<<" , thisVtx->get_x() = "<<thisVtx->get_x()<<std::endl;
	  }
        }
        else
        {
          true_daughter_vertex_x_method2 = 0;
          true_daughter_vertex_y_method2 = 0;
          true_daughter_vertex_z_method2 = 0;
        }
      }

      // project to R_EMCAL
      thisState = svtxtrack->get_state(caloRadiusEMCal);

      //int pdgid = kfp_daughter->GetPDG(); // pdgid might be reverse
      int charge = kfp_daughter->Q();
      totalcharge+=charge;

      if ((m_doLikesign==true && j==1) || (m_doLikesign==false && charge == -1))
      {
        kfp_em = kfp_daughter;
        _em_mass.push_back(kfp_daughter->GetMass());
        _em_x.push_back(kfp_daughter->GetX());
        _em_x_raw.push_back(svtxtrack->get_x());
        _em_y.push_back(kfp_daughter->GetY());
        _em_y_raw.push_back(svtxtrack->get_y());
        _em_z.push_back(kfp_daughter->GetZ());
        _em_z_raw.push_back(svtxtrack->get_z());
        _em_px.push_back(kfp_daughter->GetPx());
        _em_px_raw.push_back(svtxtrack->get_px());
        _em_py.push_back(kfp_daughter->GetPy());
        _em_py_raw.push_back(svtxtrack->get_py());
        _em_pz.push_back(kfp_daughter->GetPz());
        _em_pz_raw.push_back(svtxtrack->get_pz());
        _em_pE.push_back(kfp_daughter->GetE());
        _em_pE_unmoved.push_back(e_daughter_unmoved);
        _em_pT.push_back(kfp_daughter->GetPt());
        _em_pTErr.push_back(kfp_daughter->GetErrPt());
        _em_pT_raw.push_back(svtxtrack->get_pt());
        _em_pT_unmoved.push_back(pt_daughter_unmoved);
        _em_p.push_back(kfp_daughter->GetP());
        _em_pErr.push_back(kfp_daughter->GetErrP());
        _em_p_raw.push_back(svtxtrack->get_p());
        _em_p_unmoved.push_back(p_daughter_unmoved);
        _em_pseudorapidity.push_back(kfp_daughter->GetEta());
        _em_pseudorapidity_raw.push_back(svtxtrack->get_eta());
        _em_rapidity.push_back(kfp_daughter->GetRapidity());
        _em_theta.push_back(kfp_daughter->GetTheta());
        _em_phi.push_back(kfp_daughter->GetPhi());
        _em_phi_raw.push_back(svtxtrack->get_phi());
        _em_chi2.push_back(kfp_daughter->GetChi2());
        _em_chi2_raw.push_back(svtxtrack->get_chisq());
        _em_nDoF.push_back(kfp_daughter->GetNDF());
        _em_nDoF_raw.push_back(svtxtrack->get_ndf());
        _em_crossing.push_back(svtxtrack->get_crossing());

        if (m_doSimulation && m_doTruthMatching)
        {
          _em_has_truthmatching.push_back(isParticleValid);
          _em_true_id.push_back(true_id);
          _em_true_px.push_back(true_px);
          _em_true_py.push_back(true_py);
          _em_true_pz.push_back(true_pz);
          _em_true_vertex_x.push_back(true_daughter_vertex_x);
          _em_true_vertex_y.push_back(true_daughter_vertex_y);
          _em_true_vertex_z.push_back(true_daughter_vertex_z);
          _em_true_vertex_x_method2.push_back(true_daughter_vertex_x_method2);
          _em_true_vertex_y_method2.push_back(true_daughter_vertex_y_method2);
          _em_true_vertex_z_method2.push_back(true_daughter_vertex_z_method2);
        }

	const auto cluster_keys(get_cluster_keys(svtxtrack));
	_em_nmvtx.push_back(count_clusters<TrkrDefs::mvtxId>(cluster_keys));
	_em_nintt.push_back(count_clusters<TrkrDefs::inttId>(cluster_keys));
	_em_ntpc.push_back(count_clusters<TrkrDefs::tpcId>(cluster_keys));
	_em_ntpot.push_back(count_clusters<TrkrDefs::micromegasId>(cluster_keys));

	for(const auto& cluster_key: cluster_keys)
	{
	  trkrCluster = trkrContainer->findCluster(cluster_key);
          if(!trkrCluster)
          {
            continue;
          }
          Acts::Vector3 global(0., 0., 0.);
          global = acts_Geometry->getGlobalPosition(cluster_key, trkrCluster);
          _em_clus_ican.push_back(i);
          //_em_clus_type.push_back(TrkrDefs::getTrkrId(cluster_key));
          _em_clus_x.push_back(global[0]);
          _em_clus_y.push_back(global[1]);
          _em_clus_z.push_back(global[2]);
	}

        if(!thisState)
        {
          _em_phi_emc.push_back(NAN);
          _em_eta_emc.push_back(NAN);
          _em_px_emc.push_back(NAN);
          _em_py_emc.push_back(NAN);
          _em_pz_emc.push_back(NAN);
          _em_x_emc.push_back(NAN);
          _em_y_emc.push_back(NAN);
          _em_z_emc.push_back(NAN);
        }
        else
        {
          _em_phi_emc.push_back(thisState->get_phi());
          _em_eta_emc.push_back(thisState->get_eta());
          _em_px_emc.push_back(thisState->get_px());
          _em_py_emc.push_back(thisState->get_py());
          _em_pz_emc.push_back(thisState->get_pz());
          _em_x_emc.push_back(thisState->get_x());
          _em_y_emc.push_back(thisState->get_y());
          _em_z_emc.push_back(thisState->get_z());
        }

      }
      else if ((m_doLikesign==true && j==2) || (m_doLikesign==false && charge == 1))
      {
        kfp_ep = kfp_daughter;
        _ep_mass.push_back(kfp_daughter->GetMass());
        _ep_x.push_back(kfp_daughter->GetX());
        _ep_x_raw.push_back(svtxtrack->get_x());
        _ep_y.push_back(kfp_daughter->GetY());
        _ep_y_raw.push_back(svtxtrack->get_y());
        _ep_z.push_back(kfp_daughter->GetZ());
        _ep_z_raw.push_back(svtxtrack->get_z());
        _ep_px.push_back(kfp_daughter->GetPx());
        _ep_px_raw.push_back(svtxtrack->get_px());
        _ep_py.push_back(kfp_daughter->GetPy());
        _ep_py_raw.push_back(svtxtrack->get_py());
        _ep_pz.push_back(kfp_daughter->GetPz());
        _ep_pz_raw.push_back(svtxtrack->get_pz());
        _ep_pE.push_back(kfp_daughter->GetE());
        _ep_pE_unmoved.push_back(e_daughter_unmoved);
        _ep_pT.push_back(kfp_daughter->GetPt());
        _ep_pTErr.push_back(kfp_daughter->GetErrPt());
        _ep_pT_raw.push_back(svtxtrack->get_pt());
        _ep_pT_unmoved.push_back(pt_daughter_unmoved);
        _ep_p.push_back(kfp_daughter->GetP());
        _ep_pErr.push_back(kfp_daughter->GetErrP());
        _ep_p_raw.push_back(svtxtrack->get_p());
        _ep_p_unmoved.push_back(p_daughter_unmoved);
        _ep_pseudorapidity.push_back(kfp_daughter->GetEta());
        _ep_pseudorapidity_raw.push_back(svtxtrack->get_eta());
        _ep_rapidity.push_back(kfp_daughter->GetRapidity());
        _ep_theta.push_back(kfp_daughter->GetTheta());
        _ep_phi.push_back(kfp_daughter->GetPhi());
        _ep_phi_raw.push_back(svtxtrack->get_phi());
        _ep_chi2.push_back(kfp_daughter->GetChi2());
        _ep_chi2_raw.push_back(svtxtrack->get_chisq());
        _ep_nDoF.push_back(kfp_daughter->GetNDF());
        _ep_nDoF_raw.push_back(svtxtrack->get_ndf());
        _ep_crossing.push_back(svtxtrack->get_crossing());

        if (m_doSimulation && m_doTruthMatching)
        {
          _ep_has_truthmatching.push_back(isParticleValid);
          _ep_true_id.push_back(true_id);
          _ep_true_px.push_back(true_px);
          _ep_true_py.push_back(true_py);
          _ep_true_pz.push_back(true_pz);
          _ep_true_vertex_x.push_back(true_daughter_vertex_x);
          _ep_true_vertex_y.push_back(true_daughter_vertex_y);
          _ep_true_vertex_z.push_back(true_daughter_vertex_z);
          _ep_true_vertex_x_method2.push_back(true_daughter_vertex_x_method2);
          _ep_true_vertex_y_method2.push_back(true_daughter_vertex_y_method2);
          _ep_true_vertex_z_method2.push_back(true_daughter_vertex_z_method2);
        }

	const auto cluster_keys(get_cluster_keys(svtxtrack));
	_ep_nmvtx.push_back(count_clusters<TrkrDefs::mvtxId>(cluster_keys));
	_ep_nintt.push_back(count_clusters<TrkrDefs::inttId>(cluster_keys));
	_ep_ntpc.push_back(count_clusters<TrkrDefs::tpcId>(cluster_keys));
	_ep_ntpot.push_back(count_clusters<TrkrDefs::micromegasId>(cluster_keys));

	for(const auto& cluster_key: cluster_keys)
	{
	  trkrCluster = trkrContainer->findCluster(cluster_key);
          if(!trkrCluster)
          {
            continue;
          }
          Acts::Vector3 global(0., 0., 0.);
          global = acts_Geometry->getGlobalPosition(cluster_key, trkrCluster);
          _ep_clus_ican.push_back(i);
          //_ep_clus_type.push_back(TrkrDefs::getTrkrId(cluster_key));
          _ep_clus_x.push_back(global[0]);
          _ep_clus_y.push_back(global[1]);
          _ep_clus_z.push_back(global[2]);
	}

        if(!thisState)
        {
          _ep_phi_emc.push_back(NAN);
          _ep_eta_emc.push_back(NAN);
          _ep_px_emc.push_back(NAN);
          _ep_py_emc.push_back(NAN);
          _ep_pz_emc.push_back(NAN);
          _ep_x_emc.push_back(NAN);
          _ep_y_emc.push_back(NAN);
          _ep_z_emc.push_back(NAN);
        }
        else
        {
          _ep_phi_emc.push_back(thisState->get_phi());
          _ep_eta_emc.push_back(thisState->get_eta());
          _ep_px_emc.push_back(thisState->get_px());
          _ep_py_emc.push_back(thisState->get_py());
          _ep_pz_emc.push_back(thisState->get_pz());
          _ep_x_emc.push_back(thisState->get_x());
          _ep_y_emc.push_back(thisState->get_y());
          _ep_z_emc.push_back(thisState->get_z());
        }

      }

    }

    _epem_DCA_2d.push_back(kfp_ep->GetDistanceFromParticleXY(*kfp_em));
    _epem_DCA_3d.push_back(kfp_ep->GetDistanceFromParticle(*kfp_em));
    KFVertex motherDecayVertex;
    motherDecayVertex += *kfp_ep;
    motherDecayVertex += *kfp_em;
    _gamma_SV_chi2_per_nDoF.push_back( motherDecayVertex.GetChi2() / motherDecayVertex.GetNDF() );
    _gamma_charge.push_back(totalcharge);

  }

  RawCluster *cluster = nullptr;

  RawClusterContainer::Range begin_end_EMC = clustersEM->getClusters();
  RawClusterContainer::Iterator clusIter_EMC;

  /// Loop over the EMCal clusters
  for (clusIter_EMC = begin_end_EMC.first; clusIter_EMC != begin_end_EMC.second; ++clusIter_EMC)
  {
    cluster = clusIter_EMC->second;
    if(cluster->get_energy() < m_emcal_e_low_cut) continue;

    _emcal_e.push_back(cluster->get_energy());
    _emcal_phi.push_back(RawClusterUtility::GetAzimuthAngle(*cluster, vertex));
    _emcal_eta.push_back(RawClusterUtility::GetPseudorapidity(*cluster, vertex));
    _emcal_x.push_back(cluster->get_x());
    _emcal_y.push_back(cluster->get_y());
    _emcal_z.push_back(cluster->get_z());
  }

  if (m_doSimulation && m_doTruthMatching)
  {
    _true_numCan = m_decayMap->size();
    for (auto &iter : *m_decayMap)
    {
      Decay decay = iter.second;

      if (decay.size() % 3 != 0)
      {
        std::cout<<"Why Truth Container is not 3*n? Skip!"<<std::endl;
        return;
      }

      TLorentzVector *motherTrueLV = new TLorentzVector();
      TLorentzVector *daughterTrueLV = new TLorentzVector();
      TVector3 *mother3Vector = new TVector3();
      TVector3 *daughter3Vector = new TVector3();
      int grandmotherID = -9999;

      HepMC::GenEvent *theEvent = nullptr;
      if (m_geneventmap)
      {
        m_genevt = m_geneventmap->get(decay[0].first.first);
        assert(m_genevt);

        theEvent = m_genevt->getEvent();

        HepMC::GenParticle *motherHepMC = theEvent->barcode_to_particle(decay[0].first.second);
        assert(motherHepMC);

        motherTrueLV->SetPxPyPzE(motherHepMC->momentum().px(), motherHepMC->momentum().py(), motherHepMC->momentum().pz(), motherHepMC->momentum().e());

        // Now get the production vertex position
        HepMC::GenVertex *thisVtx = motherHepMC->production_vertex();

        if (thisVtx != nullptr) {
          mother3Vector->SetXYZ(thisVtx->point3d().x(), thisVtx->point3d().y(), thisVtx->point3d().z());
          for (auto grandmother = thisVtx->particles_in_const_begin(); grandmother != thisVtx->particles_in_const_end(); grandmother++) {
            grandmotherID = (*grandmother)->pdg_id();
            if(Verbosity() > 2) std::cout<<"HepMC grandmotherID = "<<grandmotherID<<std::endl;
          }
        }
        else {
          mother3Vector->SetXYZ(-999,-999,-999);
          if(Verbosity() > 2) std::cout<<"No grandmotherID in HepMC, because no production vertex."<<std::endl;
        }
      }

      for (unsigned int i = 1; i < decay.size(); ++i)
      {
        int pid = -999;
        if (theEvent && decay[i].first.second > -1)
        {
          HepMC::GenParticle *daughterHepMC = theEvent->barcode_to_particle(decay[i].first.second);
          daughterTrueLV->SetPxPyPzE(daughterHepMC->momentum().px(), daughterHepMC->momentum().py(), daughterHepMC->momentum().pz(), daughterHepMC->momentum().e());
          pid = daughterHepMC->pdg_id();

          // Now get the decay vertex position
          HepMC::GenVertex *thisVtx = daughterHepMC->production_vertex();

          daughter3Vector->SetXYZ(thisVtx->point3d().x(), thisVtx->point3d().y(), thisVtx->point3d().z());
        }
        else
        {
          PHG4TruthInfoContainer::ConstRange range = m_truthInfo->GetParticleRange();
          for (PHG4TruthInfoContainer::ConstIterator iter2 = range.first; iter2 != range.second; ++iter2)
          {
            PHG4Particle *daughterG4 = iter2->second;

            PHG4Particle *motherG4 = nullptr;
            PHG4Particle *grandmotherG4 = nullptr;
            if (daughterG4->get_parent_id() != 0)
            {
              motherG4 = m_truthInfo->GetParticle(daughterG4->get_parent_id());
              if (motherG4->get_parent_id() != 0)
              {
                grandmotherG4 = m_truthInfo->GetParticle(motherG4->get_parent_id());
              }
            }
            else
            {
              continue;
            }

            if (motherG4->get_pid() == decay[0].second && motherG4->get_barcode() == decay[0].first.second && daughterG4->get_pid() == decay[i].second && daughterG4->get_barcode() == decay[i].first.second)
            {
              pid = daughterG4->get_pid();

              TVector3 *motherTrue3Vector = new TVector3(motherG4->get_px(), motherG4->get_py(), motherG4->get_pz());
              motherTrueLV->SetVectM((*motherTrue3Vector), getParticleMass(decay[0].second));

              PHG4VtxPoint *thisVtx = m_truthInfo->GetVtx(motherG4->get_vtx_id());
              mother3Vector->SetXYZ(thisVtx->get_x(), thisVtx->get_y(), thisVtx->get_z());

              daughterTrueLV->SetVectM(TVector3(daughterG4->get_px(), daughterG4->get_py(), daughterG4->get_pz()), getParticleMass(decay[i].second));

              // Now get the decay vertex position
              thisVtx = m_truthInfo->GetVtx(daughterG4->get_vtx_id());
              daughter3Vector->SetXYZ(thisVtx->get_x(), thisVtx->get_y(), thisVtx->get_z());
	      if (Verbosity() > 2)
	      {
	        std::cout<<"daughterG4->get_pid() = "<<daughterG4->get_pid()<<" , daughterG4->get_barcode() = "<<daughterG4->get_barcode()<<" , thisVtx->get_x() = "<<thisVtx->get_x()<<std::endl;
	      }

              if (grandmotherG4)
              {
                grandmotherID = grandmotherG4->get_pid();
                if (Verbosity() > 2) {std::cout<<"G4 grandmotherID = "<<grandmotherID<<std::endl;}
              }
              else
              {
                if (Verbosity() > 2) {std::cout<<"No grandmotherID in G4."<<std::endl;}
              }

              delete motherTrue3Vector;
            }
          }
        }

        // e+ pdgid = -11, e- pdgid = 11
        if (pid==-11)
        {
          _true_ep_px.push_back(daughterTrueLV->Px());
          _true_ep_py.push_back(daughterTrueLV->Py());
          _true_ep_pz.push_back(daughterTrueLV->Pz());
          _true_ep_pE.push_back(daughterTrueLV->E());
          _true_ep_eta.push_back(daughterTrueLV->PseudoRapidity());
          _true_ep_phi.push_back(daughterTrueLV->Phi());
          _true_ep_x.push_back(daughter3Vector->X());
          _true_ep_y.push_back(daughter3Vector->Y());
          _true_ep_z.push_back(daughter3Vector->Z());
        }
        else if (pid==11)
        {
          _true_em_px.push_back(daughterTrueLV->Px());
          _true_em_py.push_back(daughterTrueLV->Py());
          _true_em_pz.push_back(daughterTrueLV->Pz());
          _true_em_pE.push_back(daughterTrueLV->E());
          _true_em_eta.push_back(daughterTrueLV->PseudoRapidity());
          _true_em_phi.push_back(daughterTrueLV->Phi());
          _true_em_x.push_back(daughter3Vector->X());
          _true_em_y.push_back(daughter3Vector->Y());
          _true_em_z.push_back(daughter3Vector->Z());
        }

      }

      _true_gamma_px.push_back(motherTrueLV->Px());
      _true_gamma_py.push_back(motherTrueLV->Py());
      _true_gamma_pz.push_back(motherTrueLV->Pz());
      _true_gamma_pE.push_back(motherTrueLV->E());
      _true_gamma_eta.push_back(motherTrueLV->PseudoRapidity());
      _true_gamma_phi.push_back(motherTrueLV->Phi());
      _true_gamma_x.push_back(mother3Vector->X());
      _true_gamma_y.push_back(mother3Vector->Y());
      _true_gamma_z.push_back(mother3Vector->Z());
      _true_gamma_mother_id.push_back(grandmotherID);
      _true_gamma_embedding_id.push_back(decay[0].first.first);
    }

  }

  _tree_KFP->Fill();

}

//____________________________________________________________________________..
bool TrackToCalo::checkTrack(SvtxTrack* track)
{
  if(!track)
  {
    return false;  
  }

  if(track->get_pt() < m_track_pt_low_cut)
  {
    return false;
  }

  if(track->get_quality() > m_track_quality)
  {
    return false;
  }

  const auto cluster_keys(get_cluster_keys(track));
  if (count_clusters<TrkrDefs::mvtxId>(cluster_keys) < m_nmvtx_low_cut)
  {
    return false;
  }
  if (count_clusters<TrkrDefs::inttId>(cluster_keys) < m_nintt_low_cut)
  {
    return false;
  }
  if (count_clusters<TrkrDefs::tpcId>(cluster_keys) < m_ntpc_low_cut)
  {
    return false;
  }
  if (count_clusters<TrkrDefs::micromegasId>(cluster_keys) < m_ntpot_low_cut)
  {
    return false;
  }

  return true;
}

//____________________________________________________________________________..
int TrackToCalo::End(PHCompositeNode *topNode)
{
  std::cout << topNode << std::endl;
  _outfile->cd();
  _outfile->Write();
  _outfile->Close();
  return Fun4AllReturnCodes::EVENT_OK;
}

void TrackToCalo::ResetTreeVectors()
{
  _vertex_id.clear();
  _vertex_crossing.clear();
  _vertex_ntracks.clear();
  _vertex_x.clear();
  _vertex_y.clear();
  _vertex_z.clear();
  _cluster_x.clear();
  _cluster_y.clear();
  _cluster_z.clear();
  _track_id.clear();
  _track_bc.clear();
  _track_phi.clear();
  _track_eta.clear();
  _track_pcax.clear();
  _track_pcay.clear();
  _track_pcaz.clear();
  _track_crossing.clear();
  _track_vx.clear();
  _track_vy.clear();
  _track_vz.clear();
  _track_quality.clear();
  _track_dcaxy.clear();
  _track_dcaz.clear();
  _track_nc_mvtx.clear();
  _track_nc_intt.clear();
  _track_nc_tpc.clear();
  _track_nc_tpot.clear();
  _track_ptq.clear();
  _track_px.clear();
  _track_py.clear();
  _track_pz.clear();
  _track_phi_origin.clear();
  _track_eta_origin.clear();
  _track_px_origin.clear();
  _track_py_origin.clear();
  _track_pz_origin.clear();
  _track_x_origin.clear();
  _track_y_origin.clear();
  _track_z_origin.clear();
  _track_phi_emc.clear();
  _track_eta_emc.clear();
  _track_px_emc.clear();
  _track_py_emc.clear();
  _track_pz_emc.clear();
  _track_x_emc.clear();
  _track_y_emc.clear();
  _track_z_emc.clear();
  _track_phi_ihc.clear();
  _track_eta_ihc.clear();
  _track_px_ihc.clear();
  _track_py_ihc.clear();
  _track_pz_ihc.clear();
  _track_x_ihc.clear();
  _track_y_ihc.clear();
  _track_z_ihc.clear();
  _track_phi_ohc.clear();
  _track_eta_ohc.clear();
  _track_px_ohc.clear();
  _track_py_ohc.clear();
  _track_pz_ohc.clear();
  _track_x_ohc.clear();
  _track_y_ohc.clear();
  _track_z_ohc.clear();
  _trClus_track_id.clear();
  _trClus_type.clear();
  _trClus_x.clear();
  _trClus_y.clear();
  _trClus_z.clear();
  _emcal_id.clear();
  _emcal_phi.clear();
  _emcal_eta.clear();
  _emcal_x.clear();
  _emcal_y.clear();
  _emcal_z.clear();
  _emcal_e.clear();
  _emcal_ecore.clear();
  _emcal_prob.clear();
  _emcal_chi2.clear();
  _emcal_tower_cluster_id.clear();
  _emcal_tower_e.clear();
  _emcal_tower_phi.clear();
  _emcal_tower_eta.clear();
  _emcal_tower_status.clear();
  _hcal_id.clear();
  _hcal_phi.clear();
  _hcal_eta.clear();
  _hcal_x.clear();
  _hcal_y.clear();
  _hcal_z.clear();
  _hcal_e.clear();
  _hcal_tower_cluster_id.clear();
  _hcal_tower_e.clear();
  _hcal_tower_phi.clear();
  _hcal_tower_eta.clear();
  _hcal_tower_status.clear();
  _hcal_tower_io.clear();
  _mbd_x.clear();
  _mbd_y.clear();
  _mbd_z.clear();
  _triggers.clear();
  _ntracks.clear();
}

void TrackToCalo::ResetTreeVectors_KFP()
{
  _gamma_mass.clear();
  _gamma_massErr.clear();
  _gamma_charge.clear();
  _gamma_x.clear();
  _gamma_y.clear();
  _gamma_z.clear();
  _gamma_px.clear();
  _gamma_py.clear();
  _gamma_pz.clear();
  _gamma_pE.clear();
  _gamma_pT.clear();
  _gamma_pTErr.clear();
  _gamma_p.clear();
  _gamma_pErr.clear();
  _gamma_pseudorapidity.clear();
  _gamma_rapidity.clear();
  _gamma_theta.clear();
  _gamma_phi.clear();
  _gamma_chi2.clear();
  _gamma_nDoF.clear();
  _gamma_vertex_volume.clear();
  _gamma_SV_chi2_per_nDoF.clear();
  _ep_mass.clear();
  _ep_x.clear();
  _ep_x_raw.clear();
  _ep_y.clear();
  _ep_y_raw.clear();
  _ep_z.clear();
  _ep_z_raw.clear();
  _ep_px.clear();
  _ep_px_raw.clear();
  _ep_py.clear();
  _ep_py_raw.clear();
  _ep_pz.clear();
  _ep_pz_raw.clear();
  _ep_pE.clear();
  _ep_pE_unmoved.clear();
  _ep_pT.clear();
  _ep_pTErr.clear();
  _ep_pT_raw.clear();
  _ep_pT_unmoved.clear();
  _ep_p.clear();
  _ep_pErr.clear();
  _ep_p_raw.clear();
  _ep_p_unmoved.clear();
  _ep_pseudorapidity.clear();
  _ep_pseudorapidity_raw.clear();
  _ep_rapidity.clear();
  _ep_theta.clear();
  _ep_phi.clear();
  _ep_phi_raw.clear();
  _ep_chi2.clear();
  _ep_chi2_raw.clear();
  _ep_nDoF.clear();
  _ep_nDoF_raw.clear();
  _ep_crossing.clear();
  _ep_nmvtx.clear();
  _ep_nintt.clear();
  _ep_ntpc.clear();
  _ep_ntpot.clear();
  _ep_clus_ican.clear();
  _ep_clus_x.clear();
  _ep_clus_y.clear();
  _ep_clus_z.clear();
  _ep_phi_emc.clear();
  _ep_eta_emc.clear();
  _ep_px_emc.clear();
  _ep_py_emc.clear();
  _ep_pz_emc.clear();
  _ep_x_emc.clear();
  _ep_y_emc.clear();
  _ep_z_emc.clear();
  _em_mass.clear();
  _em_x.clear();
  _em_x_raw.clear();
  _em_y.clear();
  _em_y_raw.clear();
  _em_z.clear();
  _em_z_raw.clear();
  _em_px.clear();
  _em_px_raw.clear();
  _em_py.clear();
  _em_py_raw.clear();
  _em_pz.clear();
  _em_pz_raw.clear();
  _em_pE.clear();
  _em_pE_unmoved.clear();
  _em_pT.clear();
  _em_pTErr.clear();
  _em_pT_raw.clear();
  _em_pT_unmoved.clear();
  _em_p.clear();
  _em_pErr.clear();
  _em_p_raw.clear();
  _em_p_unmoved.clear();
  _em_pseudorapidity.clear();
  _em_pseudorapidity_raw.clear();
  _em_rapidity.clear();
  _em_theta.clear();
  _em_phi.clear();
  _em_phi_raw.clear();
  _em_chi2.clear();
  _em_chi2_raw.clear();
  _em_nDoF.clear();
  _em_nDoF_raw.clear();
  _em_crossing.clear();
  _em_nmvtx.clear();
  _em_nintt.clear();
  _em_ntpc.clear();
  _em_ntpot.clear();
  _em_clus_ican.clear();
  _em_clus_x.clear();
  _em_clus_y.clear();
  _em_clus_z.clear();
  _em_phi_emc.clear();
  _em_eta_emc.clear();
  _em_px_emc.clear();
  _em_py_emc.clear();
  _em_pz_emc.clear();
  _em_x_emc.clear();
  _em_y_emc.clear();
  _em_z_emc.clear();
  _emcal_phi.clear();
  _emcal_eta.clear();
  _emcal_x.clear();
  _emcal_y.clear();
  _emcal_z.clear();
  _emcal_e.clear();
  _epem_DCA_2d.clear();
  _epem_DCA_3d.clear();

  _em_has_truthmatching.clear();
  _em_true_id.clear();
  _em_true_px.clear();
  _em_true_py.clear();
  _em_true_pz.clear();
  _em_true_vertex_x.clear();
  _em_true_vertex_y.clear();
  _em_true_vertex_z.clear();
  _em_true_vertex_x_method2.clear();
  _em_true_vertex_y_method2.clear();
  _em_true_vertex_z_method2.clear();

  _ep_has_truthmatching.clear();
  _ep_true_id.clear();
  _ep_true_px.clear();
  _ep_true_py.clear();
  _ep_true_pz.clear();
  _ep_true_vertex_x.clear();
  _ep_true_vertex_y.clear();
  _ep_true_vertex_z.clear();
  _ep_true_vertex_x_method2.clear();
  _ep_true_vertex_y_method2.clear();
  _ep_true_vertex_z_method2.clear();

  _true_gamma_phi.clear();
  _true_gamma_eta.clear();
  _true_gamma_px.clear();
  _true_gamma_py.clear();
  _true_gamma_pz.clear();
  _true_gamma_pE.clear();
  _true_gamma_x.clear();
  _true_gamma_y.clear();
  _true_gamma_z.clear();
  _true_gamma_mother_id.clear();
  _true_gamma_embedding_id.clear();

  _true_ep_phi.clear();
  _true_ep_eta.clear();
  _true_ep_px.clear();
  _true_ep_py.clear();
  _true_ep_pz.clear();
  _true_ep_pE.clear();
  _true_ep_x.clear();
  _true_ep_y.clear();
  _true_ep_z.clear();

  _true_em_phi.clear();
  _true_em_eta.clear();
  _true_em_px.clear();
  _true_em_py.clear();
  _true_em_pz.clear();
  _true_em_pE.clear();
  _true_em_x.clear();
  _true_em_y.clear();
  _true_em_z.clear();

}

void TrackToCalo::resetCaloRadius()
{
  if (!EMCalGeo || !IHCalGeo || !OHCalGeo)
  {
    std::cout << PHWHERE << "missing node trees, can't continue with track projection"
              << std::endl;
    return;
  }

  if (m_use_emcal_radius)
  {
    caloRadiusEMCal = m_emcal_radius_user;
  }
  else
  {
    caloRadiusEMCal = EMCalGeo->get_radius();
  }
  if (m_use_ihcal_radius)
  {
    caloRadiusIHCal = m_ihcal_radius_user;
  }
  else
  {
    caloRadiusIHCal = IHCalGeo->get_radius();
  }
  if (m_use_ohcal_radius)
  {
    caloRadiusOHCal = m_ohcal_radius_user;
  }
  else
  {
    caloRadiusOHCal = OHCalGeo->get_radius();
  }
}

PHG4Particle *TrackToCalo::getTruthTrack(SvtxTrack *thisTrack)
{
  /*
   * There are two methods for getting the truth rack from the reco track
   * 1. (recommended) Use the reco -> truth tables (requires SvtxPHG4ParticleMap). Introduced Summer of 2022
   * 2. Get truth track via nClusters. Older method and will work with older DSTs
   */

  PHG4Particle *particle = nullptr;

  if (dst_reco_truth_map)
  {
    std::map<float, std::set<int>> truth_set = dst_reco_truth_map->get(thisTrack->get_id());
    const auto &best_weight = truth_set.rbegin();
    int best_truth_id = *best_weight->second.rbegin();
    particle = m_truthInfo->GetParticle(best_truth_id);
  }
  else
  {
    std::cout << __FILE__ << ": SvtxPHG4ParticleMap not found, reverting to max_truth_particle_by_nclusters()" << std::endl;
    particle = trackeval->max_truth_particle_by_nclusters(thisTrack);
  }
  return particle;
}
