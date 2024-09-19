// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/*!
 *  \file               TrackToCalo.h
 *  \brief              Track To Calo, output root file
 *  \author Antonio Silva <antonio.silva@cern.ch>, Xudong Yu <xyu3@bnl.gov>
 */

#ifndef TRACKTOCALO_H
#define TRACKTOCALO_H

#include <fun4all/SubsysReco.h>

#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/SvtxVertexMap.h>
#include <ffarawobjects/Gl1Packet.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/TrackSeed.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrCluster.h>
#include <tpc/LaserEventInfo.h>
#include <calobase/RawTowerGeomContainer.h>

#include <kfparticle_sphenix/KFParticle_sPHENIX.h>
#include <kfparticle_sphenix/KFParticle_DST.h>
#include <kfparticle_sphenix/KFParticle_Container.h>
#include <kfparticle_sphenix/KFParticle_Tools.h>
#include <kfparticle_sphenix/KFParticle_truthAndDetTools.h>
#include <kfparticle_sphenix/KFParticle_nTuple.h>

#include <string>
#include <vector>

class PHCompositeNode;
class TH1;
class TH2;
class TFile;
class TTree;
class KFParticle_sPHENIX;

class TrackToCalo : public SubsysReco
{
 public:

  TrackToCalo(const std::string &name = "TrackToCalo", const std::string &file = "output.root");

  ~TrackToCalo() override;

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  void ResetTreeVectors();
  void ResetTreeVectors_KFP();

  void fillTree();
  void fillTree_TrackOnly();
  void fillTree_CaloOnly();
  void fillTree_KFP();

  void createBranches();
  void createBranches_KFP();

  void EMcalRadiusUser(bool use) {m_use_emcal_radius = use;}
  void IHcalRadiusUser(bool use) {m_use_ihcal_radius = use;}
  void OHcalRadiusUser(bool use) {m_use_ohcal_radius = use;}
  void setEMcalRadius(float r) {m_emcal_radius_user = r;}
  void setIHcalRadius(float r) {m_ihcal_radius_user = r;}
  void setOHcalRadius(float r) {m_ohcal_radius_user = r;}

  void setRawClusContEMName(std::string name) {m_RawClusCont_EM_name = name;}
  void setRawClusContHADName(std::string name) {m_RawClusCont_HAD_name = name;}
  void setKFPContName(std::string name) {m_KFPCont_name = name;}
  void setKFPtrackMapName(std::string name) {m_KFPtrackMap_name = name;}

  void resetCaloRadius();

  void setTrackPtLowCut(float pt) {m_track_pt_low_cut = pt;}
  void setEmcalELowCut(float e) {m_emcal_e_low_cut = e;}

  void doTrkrCaloMatching() {m_doTrkrCaloMatching = true;}
  void doTrkrCaloMatching_KFP() {m_doTrkrCaloMatching_KFP = true;}

  void anaTrkrInfo() {m_doTrackOnly = true;}
  void anaCaloInfo() {m_doCaloOnly = true;}

  void setRejectLaserEvent(bool value) {m_rejectLaserEvent = value;}

 private:
   int cnt = 0;
   bool m_use_emcal_radius = false;
   bool m_use_ihcal_radius = false;
   bool m_use_ohcal_radius = false;
   float m_emcal_radius_user = 93.5;
   float m_ihcal_radius_user = 117;
   float m_ohcal_radius_user = 177.423;
   std::string _outfilename;
   TFile *_outfile = nullptr;
   TTree *_tree = nullptr;
   TTree *_tree_KFP = nullptr;

   std::string m_RawClusCont_EM_name = "TOPOCLUSTER_EMCAL";
   std::string m_RawClusCont_HAD_name = "TOPOCLUSTER_HCAL";

   std::string m_KFPCont_name = "KFParticle_Container";
   std::string m_KFPtrackMap_name = "SvtxTrackMap";

   int _runNumber;
   int _eventNumber;
   std::vector<int> _vertex_id;
   std::vector<int> _vertex_crossing;
   std::vector<int> _vertex_ntracks;
   std::vector<float> _vertex_x;
   std::vector<float> _vertex_y;
   std::vector<float> _vertex_z;
   std::vector<float> _cluster_x;
   std::vector<float> _cluster_y;
   std::vector<float> _cluster_z;
   std::vector<int> _track_id;
   std::vector<int> _track_bc;
   std::vector<float> _track_phi;
   std::vector<float> _track_eta;
   std::vector<float> _track_pcax;
   std::vector<float> _track_pcay;
   std::vector<float> _track_pcaz;
   std::vector<float> _track_crossing;
   std::vector<float> _track_vx;
   std::vector<float> _track_vy;
   std::vector<float> _track_vz;
   std::vector<float> _track_quality;
   std::vector<float> _track_dcaxy;
   std::vector<float> _track_dcaz;
   std::vector<int> _track_nc_mvtx;
   std::vector<int> _track_nc_intt;
   std::vector<int> _track_nc_tpc;
   std::vector<float> _track_ptq;
   std::vector<float> _track_px;
   std::vector<float> _track_py;
   std::vector<float> _track_pz;
   std::vector<float> _track_phi_origin;
   std::vector<float> _track_eta_origin;
   std::vector<float> _track_px_origin;
   std::vector<float> _track_py_origin;
   std::vector<float> _track_pz_origin;
   std::vector<float> _track_x_origin;
   std::vector<float> _track_y_origin;
   std::vector<float> _track_z_origin;
   std::vector<float> _track_phi_emc;
   std::vector<float> _track_eta_emc;
   std::vector<float> _track_px_emc;
   std::vector<float> _track_py_emc;
   std::vector<float> _track_pz_emc;
   std::vector<float> _track_x_emc;
   std::vector<float> _track_y_emc;
   std::vector<float> _track_z_emc;
   std::vector<float> _track_phi_ihc;
   std::vector<float> _track_eta_ihc;
   std::vector<float> _track_px_ihc;
   std::vector<float> _track_py_ihc;
   std::vector<float> _track_pz_ihc;
   std::vector<float> _track_x_ihc;
   std::vector<float> _track_y_ihc;
   std::vector<float> _track_z_ihc;
   std::vector<float> _track_phi_ohc;
   std::vector<float> _track_eta_ohc;
   std::vector<float> _track_px_ohc;
   std::vector<float> _track_py_ohc;
   std::vector<float> _track_pz_ohc;
   std::vector<float> _track_x_ohc;
   std::vector<float> _track_y_ohc;
   std::vector<float> _track_z_ohc;

   std::vector<int> _trClus_track_id;
   std::vector<int> _trClus_type;
   std::vector<float> _trClus_x;
   std::vector<float> _trClus_y;
   std::vector<float> _trClus_z;

   std::vector<int> _emcal_id;
   std::vector<float> _emcal_phi;
   std::vector<float> _emcal_eta;
   std::vector<float> _emcal_x;
   std::vector<float> _emcal_y;
   std::vector<float> _emcal_z;
   std::vector<float> _emcal_e;
   std::vector<float> _emcal_ecore;
   std::vector<float> _emcal_chi2;
   std::vector<float> _emcal_prob;
   std::vector<int> _emcal_tower_cluster_id;
   std::vector<float> _emcal_tower_e;
   std::vector<float> _emcal_tower_phi;
   std::vector<float> _emcal_tower_eta;
   std::vector<int> _emcal_tower_status;

   std::vector<int> _hcal_id;
   std::vector<float> _hcal_phi;
   std::vector<float> _hcal_eta;
   std::vector<float> _hcal_x;
   std::vector<float> _hcal_y;
   std::vector<float> _hcal_z;
   std::vector<float> _hcal_e;
   std::vector<int> _hcal_tower_cluster_id;
   std::vector<float> _hcal_tower_e;
   std::vector<float> _hcal_tower_phi;
   std::vector<float> _hcal_tower_eta;
   std::vector<int> _hcal_tower_status;
   std::vector<int> _hcal_tower_io;

   std::vector<float> _mbd_x;
   std::vector<float> _mbd_y;
   std::vector<float> _mbd_z;

   std::vector<int> _triggers;

   std::vector<int> _ntracks;

   int _numCan;

   std::vector<float> _gamma_mass;
   std::vector<float> _gamma_massErr;
   std::vector<float> _gamma_x;
   std::vector<float> _gamma_y;
   std::vector<float> _gamma_z;
   std::vector<float> _gamma_px;
   std::vector<float> _gamma_py;
   std::vector<float> _gamma_pz;
   std::vector<float> _gamma_pE;
   std::vector<float> _gamma_pT;
   std::vector<float> _gamma_pTErr;
   std::vector<float> _gamma_p;
   std::vector<float> _gamma_pErr;
   std::vector<float> _gamma_pseudorapidity;
   std::vector<float> _gamma_rapidity;
   std::vector<float> _gamma_theta;
   std::vector<float> _gamma_phi;
   std::vector<float> _gamma_chi2;
   std::vector<float> _gamma_nDoF;
   std::vector<float> _gamma_vertex_volume;
   std::vector<float> _gamma_SV_chi2_per_nDoF;

   std::vector<float> _ep_mass;
   std::vector<float> _ep_x;
   std::vector<float> _ep_y;
   std::vector<float> _ep_z;
   std::vector<float> _ep_px;
   std::vector<float> _ep_py;
   std::vector<float> _ep_pz;
   std::vector<float> _ep_pE;
   std::vector<float> _ep_pE_unmoved;
   std::vector<float> _ep_pT;
   std::vector<float> _ep_pTErr;
   std::vector<float> _ep_pT_raw;
   std::vector<float> _ep_pT_unmoved;
   std::vector<float> _ep_p;
   std::vector<float> _ep_pErr;
   std::vector<float> _ep_p_raw;
   std::vector<float> _ep_p_unmoved;
   std::vector<float> _ep_pseudorapidity;
   std::vector<float> _ep_rapidity;
   std::vector<float> _ep_theta;
   std::vector<float> _ep_phi;
   std::vector<float> _ep_chi2;
   std::vector<float> _ep_nDoF;
   std::vector<float> _ep_crossing;
   std::vector<int> _ep_clus_ican;
   //std::vector<int> _ep_clus_type;
   std::vector<float> _ep_clus_x;
   std::vector<float> _ep_clus_y;
   std::vector<float> _ep_clus_z;

   std::vector<float> _em_mass;
   std::vector<float> _em_x;
   std::vector<float> _em_y;
   std::vector<float> _em_z;
   std::vector<float> _em_px;
   std::vector<float> _em_py;
   std::vector<float> _em_pz;
   std::vector<float> _em_pE;
   std::vector<float> _em_pE_unmoved;
   std::vector<float> _em_pT;
   std::vector<float> _em_pTErr;
   std::vector<float> _em_pT_raw;
   std::vector<float> _em_pT_unmoved;
   std::vector<float> _em_p;
   std::vector<float> _em_pErr;
   std::vector<float> _em_p_raw;
   std::vector<float> _em_p_unmoved;
   std::vector<float> _em_pseudorapidity;
   std::vector<float> _em_rapidity;
   std::vector<float> _em_theta;
   std::vector<float> _em_phi;
   std::vector<float> _em_chi2;
   std::vector<float> _em_nDoF;
   std::vector<float> _em_crossing;
   std::vector<int> _em_clus_ican;
   //std::vector<int> _em_clus_type;
   std::vector<float> _em_clus_x;
   std::vector<float> _em_clus_y;
   std::vector<float> _em_clus_z;

   std::vector<float> _ep_phi_emc;
   std::vector<float> _ep_eta_emc;
   std::vector<float> _ep_px_emc;
   std::vector<float> _ep_py_emc;
   std::vector<float> _ep_pz_emc;
   std::vector<float> _ep_x_emc;
   std::vector<float> _ep_y_emc;
   std::vector<float> _ep_z_emc;

   std::vector<float> _em_phi_emc;
   std::vector<float> _em_eta_emc;
   std::vector<float> _em_px_emc;
   std::vector<float> _em_py_emc;
   std::vector<float> _em_pz_emc;
   std::vector<float> _em_x_emc;
   std::vector<float> _em_y_emc;
   std::vector<float> _em_z_emc;

   std::vector<float> _epem_DCA_2d;
   std::vector<float> _epem_DCA_3d;

   GlobalVertexMap *vertexmap = nullptr;
   SvtxVertexMap *vertexMap = nullptr;
   Gl1Packet *gl1Packet = nullptr;
   SvtxTrackMap *trackMap = nullptr;
   SvtxTrackMap *KFP_trackMap = nullptr;
   KFParticle_Container *KFP_Container = nullptr;
   ActsGeometry *acts_Geometry = nullptr;
   RawClusterContainer *clustersEM = nullptr;
   RawClusterContainer *clustersHAD = nullptr;
   RawClusterContainer *EMCAL_RawClusters = nullptr;
   TowerInfoContainer *EMCAL_Container = nullptr;
   TowerInfoContainer *IHCAL_Container = nullptr;
   TowerInfoContainer *OHCAL_Container = nullptr;
   TrkrHitSetContainer *trkrHitSet = nullptr;
   TrkrClusterContainer *trkrContainer = nullptr;
   RawTowerGeomContainer *EMCalGeo = nullptr;
   RawTowerGeomContainer *IHCalGeo = nullptr;
   RawTowerGeomContainer *OHCalGeo = nullptr;
   LaserEventInfo* laserEventInfo = nullptr;

   SvtxTrackState *thisState = nullptr;
   SvtxTrack *track = nullptr;
   TrackSeed *seed = nullptr;
   TrackSeed *tpc_seed = nullptr;
   TrkrCluster *trkrCluster = nullptr;

   KFParticle* kfp_mother = nullptr;
   KFParticle* kfp_daughter = nullptr;
   KFParticle* kfp_ep = nullptr;
   KFParticle* kfp_em = nullptr;

   float m_track_pt_low_cut = 0.5;
   float m_emcal_e_low_cut = 0.2;
   float m_vx, m_vy, m_vz;

   double caloRadiusEMCal;
   double caloRadiusIHCal;
   double caloRadiusOHCal;

   bool m_doTrkrCaloMatching = false;
   bool m_doTrkrCaloMatching_KFP = false;
   bool m_doTrackOnly = false;
   bool m_doCaloOnly = false;

   bool m_rejectLaserEvent = true;
};

#endif // TRACKTOCALO_H
